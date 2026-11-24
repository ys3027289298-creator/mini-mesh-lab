#define _POSIX_C_SOURCE 200809L

/* interrupt: cancellable batched mesh pipeline (weld -> normals, bbox tracked
 * alongside). Cancellation is cooperative: the cancel flag is observed ONLY at
 * batch boundaries, and every batch commits atomically, so an interrupted job
 * always exposes a self-consistent state and can be resumed to the same result
 * as an uninterrupted run.
 *
 * Build:
 *   make            production binary, runs the random-cancel test suite
 *   make fault      three fault-injection binaries, each breaking one step's
 *                   atomicity; the consistency validator must name that step
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <stdarg.h>

#ifndef FAULT_INJECT
#define FAULT_INJECT 0          /* 0 clean, 1 weld compaction, 2 normal publish, 3 bbox order */
#endif

#define MAXV 1024
#define MAXF 2048
#define VBATCH 8                /* vertices welded per batch */
#define FBATCH 12               /* faces remapped / normals accumulated per batch */
#define EPS 1e-3f

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { V3 min, max; int valid; } Box;

static V3 sub3(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 cross3(V3 p, V3 q) { V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r; }
static float dot3(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float len3(V3 p) { return sqrtf(dot3(p,p)); }
static float dist3(V3 p, V3 q) { return len3(sub3(p, q)); }

static void box_empty(Box* b) {
    b->min.x = b->min.y = b->min.z = 1e30f;
    b->max.x = b->max.y = b->max.z = -1e30f;
    b->valid = 0;
}
static void box_grow(Box* b, V3 p) {
    if (b->min.x > p.x) b->min.x = p.x;
    if (b->min.y > p.y) b->min.y = p.y;
    if (b->min.z > p.z) b->min.z = p.z;
    if (b->max.x < p.x) b->max.x = p.x;
    if (b->max.y < p.y) b->max.y = p.y;
    if (b->max.z < p.z) b->max.z = p.z;
    b->valid = 1;
}

/* ---------------- input mesh ---------------- */

typedef struct { V3 v[MAXV]; int nv; Tri f[MAXF]; int nf; } MeshIn;

static int in_add_v(MeshIn* m, float x, float y, float z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv].x = x; m->v[m->nv].y = y; m->v[m->nv].z = z;
    return m->nv++;
}
static int in_add_f(MeshIn* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* ---------------- cancellable job ---------------- */

enum { ST_WELD = 0, ST_NORMALS = 1, ST_DONE = 2 };

typedef struct {
    const MeshIn* in;
    int stage;                  /* ST_WELD -> ST_NORMALS -> ST_DONE */
    int fault_armed;            /* fault demo: cancel at checkpoint fault_cp */
    int fault_cp;
    int sig_arm;                /* test scaffold: raise(SIGINT) at this cp, -1 = off */

    /* weld stage cursors: everything below them is committed and consistent */
    int vcur;                   /* input vertices welded */
    int fcur;                   /* input faces ingested */
    int ncommitted;             /* non-degenerate remapped faces published */
    int nu;                     /* welded (unique) vertices so far */
    V3 uniq[MAXV];
    int map[MAXV];              /* input index -> welded index, -1 until welded */
    Tri f[MAXF];                /* committed remapped faces only */
    Box box;

    /* normal stage: acc[] is PRIVATE until the whole finalize publishes nrm[] */
    int acur;                   /* faces accumulated into acc[] */
    V3 acc[MAXV];
    V3 nrm[MAXV];
    int normals_published;
    int bad;                    /* zero-length accumulator count, set at publish */
    int vnorm;                  /* fault-2 only: per-vertex publish cursor */

    volatile sig_atomic_t cancel_request;
    int checkpoint_count;       /* batch boundaries observed */
} Job;

static volatile sig_atomic_t sig_flag;
static void on_signal(int s) { (void)s; sig_flag = 1; }

/* Install a PERSISTENT handler. Plain signal() is unsafe here: under
 * -std=c11 / POSIX feature macros glibc gives it System V one-shot semantics,
 * i.e. the disposition resets to SIG_DFL after the first delivery and the
 * second cancellation would kill the process. sigaction has no such quirk. */
static void install_handler(void) {
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = on_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;                 /* no SA_RESETHAND: stays installed */
    sigaction(SIGINT, &sa, NULL);
}

static void job_init(Job* j, const MeshIn* in) {
    memset(j, 0, sizeof(*j));
    j->in = in;
    j->stage = ST_WELD;
    j->fault_cp = -1;
    j->sig_arm = -1;
    for (int i = 0; i < MAXV; i++) j->map[i] = -1;
    box_empty(&j->box);
}

/* The single cancellation point. Called exactly once per finished batch, never
 * in the middle of one. Consumes the request so a later resume starts clean. */
static int checkpoint(Job* j) {
    int ord = j->checkpoint_count++;
    if (j->sig_arm == ord) { j->sig_arm = -1; raise(SIGINT); }
    if (j->cancel_request || sig_flag) {
        j->cancel_request = 0;
        sig_flag = 0;
        return -1;
    }
    if (j->fault_armed && ord == j->fault_cp) {
        j->fault_armed = 0;
        return -1;
    }
    return 0;
}

/* ---- stage 1a: weld one batch of input vertices (bbox grown in same step) ---- */

static int weld_lookup(Job* j, V3 p) {
    for (int k = 0; k < j->nu; k++) {
        if (dist3(p, j->uniq[k]) <= EPS) return k;
    }
    int k = j->nu++;
    j->uniq[k] = p;
    return k;
}

static int step_vertices(Job* j) {
    int end = j->vcur + VBATCH;
    if (end > j->in->nv) end = j->in->nv;
    while (j->vcur < end) {
        int i = j->vcur;
#if FAULT_INJECT == 3
        /* BUG (bbox step): grow the box from a shuffled vertex order, so at any
         * boundary the box does not describe the prefix that was welded. */
        box_grow(&j->box, j->in->v[j->in->nv - 1 - i]);
#else
        box_grow(&j->box, j->in->v[i]);
#endif
        j->map[i] = weld_lookup(j, j->in->v[i]);
        j->vcur++;
    }
    return checkpoint(j);
}

/* ---- stage 1b: remap one batch of faces, publish survivors atomically ---- */

static int step_faces(Job* j) {
    int end = j->fcur + FBATCH;
    if (end > j->in->nf) end = j->in->nf;
#if FAULT_INJECT == 1
    /* BUG (weld step): remap faces in place, over the live array, and count
     * every input face as committed. Degenerate faces and half-remapped state
     * become visible at the batch boundary -- compaction is not atomic. */
    while (j->fcur < end) {
        int i = j->fcur;
        Tri t = j->in->f[i];
        j->f[i].a = j->map[t.a]; j->f[i].b = j->map[t.b]; j->f[i].c = j->map[t.c];
        j->fcur++;
    }
    j->ncommitted = j->fcur;
#else
    Tri staged[FBATCH];
    int keep = 0;
    while (j->fcur < end) {
        int i = j->fcur;
        Tri t = j->in->f[i];
        Tri r = { j->map[t.a], j->map[t.b], j->map[t.c] };
        if (r.a != r.b && r.b != r.c && r.a != r.c) staged[keep++] = r;
        j->fcur++;
    }
    for (int k = 0; k < keep; k++) j->f[j->ncommitted + k] = staged[k];  /* atomic publish */
    j->ncommitted += keep;
#endif
    return checkpoint(j);
}

/* ---- stage 2a: accumulate unnormalized face normals, one batch at a time ---- */

static int step_accum(Job* j) {
    int end = j->acur + FBATCH;
    if (end > j->ncommitted) end = j->ncommitted;
    while (j->acur < end) {
        Tri t = j->f[j->acur];
        V3 n = cross3(sub3(j->uniq[t.b], j->uniq[t.a]),
                      sub3(j->uniq[t.c], j->uniq[t.a]));
        j->acc[t.a].x += n.x; j->acc[t.a].y += n.y; j->acc[t.a].z += n.z;
        j->acc[t.b].x += n.x; j->acc[t.b].y += n.y; j->acc[t.b].z += n.z;
        j->acc[t.c].x += n.x; j->acc[t.c].y += n.y; j->acc[t.c].z += n.z;
        j->acur++;
    }
    /* acc[] is private scratch: a cancel here reveals no half-normal state,
     * because nrm[] stays untouched and normals_published stays 0. */
    return checkpoint(j);
}

/* ---- stage 2b: normalize every vertex and publish in one non-interruptible step ---- */

static int finalize_normals(Job* j) {
#if FAULT_INJECT == 2
    /* BUG (normal step): publish per vertex instead of one atomic finalize.
     * A cancel between vertices leaves some vertices with unit normals and the
     * rest at zero while the published gate is already on. */
    j->normals_published = 1;
    while (j->vnorm < j->nu) {
        int i = j->vnorm;
        float l = len3(j->acc[i]);
        if (l >= 1e-12f) { j->nrm[i].x = j->acc[i].x/l; j->nrm[i].y = j->acc[i].y/l; j->nrm[i].z = j->acc[i].z/l; }
        j->vnorm++;
        if (j->vnorm % VBATCH == 0) {
            int r = checkpoint(j);
            if (r < 0) return r;
        }
    }
    j->bad = 0;
    for (int i = 0; i < j->nu; i++) if (len3(j->acc[i]) < 1e-12f) j->bad++;
    j->stage = ST_DONE;
    return checkpoint(j);
#else
    int bad = 0;
    for (int i = 0; i < j->nu; i++) {
        float l = len3(j->acc[i]);
        if (l < 1e-12f) { bad++; continue; }
        j->nrm[i].x = j->acc[i].x/l; j->nrm[i].y = j->acc[i].y/l; j->nrm[i].z = j->acc[i].z/l;
    }
    j->bad = bad;
    j->normals_published = 1;      /* gate flips once, after every vertex is done */
    j->stage = ST_DONE;
    return checkpoint(j);
#endif
}

/* Returns 0 while active, 1 when fully done, -1 when a batch boundary observed
 * the cancel flag. No state changes survive a return of -1 beyond the last
 * fully committed batch. */
int job_step(Job* j) {
    if (j->stage == ST_DONE) return 1;
    int r;
    if (j->stage == ST_WELD) {
        if (j->vcur < j->in->nv) r = step_vertices(j);
        else if (j->fcur < j->in->nf) r = step_faces(j);
        else { j->stage = ST_NORMALS; return 0; }
        if (r < 0) return -1;
        if (j->stage == ST_WELD) return 0;
    }
    if (j->acur < j->ncommitted) {
        r = step_accum(j);
        if (r < 0) return -1;
        return 0;
    }
    r = finalize_normals(j);
    if (r < 0) return -1;
    return 1;
}

/* ---------------- independent reference implementation ---------------- */

typedef struct {
    V3 v[MAXV]; int nv;
    Tri f[MAXF]; int nf;
    int map[MAXV];
    V3 nrm[MAXV]; int bad;
    Box box;
} Reference;

static void reference_run(const MeshIn* in, Reference* ref) {
    memset(ref, 0, sizeof(*ref));
    box_empty(&ref->box);
    for (int i = 0; i < in->nv; i++) {
        box_grow(&ref->box, in->v[i]);
        int found = -1;
        for (int k = 0; k < ref->nv; k++) {
            if (dist3(in->v[i], ref->v[k]) <= EPS) { found = k; break; }
        }
        if (found < 0) { found = ref->nv++; ref->v[found] = in->v[i]; }
        ref->map[i] = found;
    }
    for (int i = 0; i < in->nf; i++) {
        Tri t = in->f[i];
        Tri r = { ref->map[t.a], ref->map[t.b], ref->map[t.c] };
        if (r.a == r.b || r.b == r.c || r.a == r.c) continue;
        ref->f[ref->nf++] = r;
    }
    V3 acc[MAXV];
    memset(acc, 0, sizeof(acc));
    for (int i = 0; i < ref->nf; i++) {
        Tri t = ref->f[i];
        V3 n = cross3(sub3(ref->v[t.b], ref->v[t.a]), sub3(ref->v[t.c], ref->v[t.a]));
        acc[t.a].x += n.x; acc[t.a].y += n.y; acc[t.a].z += n.z;
        acc[t.b].x += n.x; acc[t.b].y += n.y; acc[t.b].z += n.z;
        acc[t.c].x += n.x; acc[t.c].y += n.y; acc[t.c].z += n.z;
    }
    for (int i = 0; i < ref->nv; i++) {
        float l = len3(acc[i]);
        if (l < 1e-12f) { ref->bad++; continue; }
        ref->nrm[i].x = acc[i].x/l; ref->nrm[i].y = acc[i].y/l; ref->nrm[i].z = acc[i].z/l;
    }
}

#if FAULT_INJECT == 0
static int tri_cmp(const void* pa, const void* pb) {
    const Tri* a = pa; const Tri* b = pb;
    if (a->a != b->a) return a->a - b->a;
    if (a->b != b->b) return a->b - b->b;
    return a->c - b->c;
}
#endif

/* ---------------- consistency validator ----------------
 * Every invariant below must hold at EVERY observable state, including the one
 * left behind by a cancellation. A failure names the processing step whose
 * batch boundary was not actually atomic. */

static int vfail(char* detail, const char* step, const char* fmt, ...) {
    char body[200];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(body, sizeof(body), fmt, ap);
    va_end(ap);
    snprintf(detail, 256, "%s: %s", step, body);
    return 0;
}

static int box_equal(const Box* a, const Box* b) {
    if (a->valid != b->valid) return 0;
    if (!a->valid) return 1;
    return fabsf(a->min.x-b->min.x) < 1e-6f && fabsf(a->min.y-b->min.y) < 1e-6f &&
           fabsf(a->min.z-b->min.z) < 1e-6f && fabsf(a->max.x-b->max.x) < 1e-6f &&
           fabsf(a->max.y-b->max.y) < 1e-6f && fabsf(a->max.z-b->max.z) < 1e-6f;
}

static int validate(const Job* j, const Reference* ref, char* detail) {
#if FAULT_INJECT != 0
    (void)ref;
#endif
    const MeshIn* in = j->in;
    detail[0] = 0;

    if (j->vcur < 0 || j->vcur > in->nv || j->fcur < 0 || j->fcur > in->nf ||
        j->nu < 0 || j->nu > j->vcur || j->ncommitted < 0 || j->ncommitted > j->fcur ||
        j->acur < 0 || j->acur > j->ncommitted)
        return vfail(detail, "cursor", "cursors out of range");

    /* weld step: map prefix must be in range and point at an equal-position vertex */
    for (int i = 0; i < j->vcur; i++) {
        int m = j->map[i];
        if (m < 0 || m >= j->nu)
            return vfail(detail, "weld.map", "vertex %d mapped outside welded set", i);
        if (dist3(in->v[i], j->uniq[m]) > EPS)
            return vfail(detail, "weld.map", "vertex %d maps to a different position", i);
    }
    /* weld step: published welded vertices must be pairwise distinct */
    for (int a = 0; a < j->nu; a++)
        for (int b = a + 1; b < j->nu; b++)
            if (dist3(j->uniq[a], j->uniq[b]) <= EPS)
                return vfail(detail, "weld.uniq", "duplicate welded vertices %d/%d", a, b);

    /* weld step: every committed face must be in welded range and non-degenerate */
    for (int i = 0; i < j->ncommitted; i++) {
        Tri t = j->f[i];
        if (t.a < 0 || t.b < 0 || t.c < 0 || t.a >= j->nu || t.b >= j->nu || t.c >= j->nu)
            return vfail(detail, "weld.face", "committed face %d indexes outside welded set", i);
        if (t.a == t.b || t.b == t.c || t.a == t.c)
            return vfail(detail, "weld.face", "committed face %d is degenerate (remap/compaction not atomic)", i);
    }

    /* bbox step: box must equal the AABB of exactly the welded prefix */
    Box expect;
    box_empty(&expect);
    int prefix = j->stage == ST_DONE ? in->nv : j->vcur;
    for (int i = 0; i < prefix; i++) box_grow(&expect, in->v[i]);
    if (!box_equal(&j->box, &expect))
        return vfail(detail, "bbox", "box does not match the %d processed vertices", prefix);

#if FAULT_INJECT == 0
    if (j->stage >= ST_NORMALS) {
        /* weld fully committed: must match the reference weld exactly */
        if (j->nu != ref->nv)
            return vfail(detail, "weld.uniq", "welded count %d != reference %d", j->nu, ref->nv);
        for (int i = 0; i < in->nv; i++)
            if (j->map[i] != ref->map[i])
                return vfail(detail, "weld.map", "map[%d]=%d != reference %d", i, j->map[i], ref->map[i]);
        if (j->ncommitted != ref->nf)
            return vfail(detail, "weld.face", "face count %d != reference %d", j->ncommitted, ref->nf);
        Tri a[MAXF], b[MAXF];
        memcpy(a, j->f, sizeof(Tri) * (size_t)j->ncommitted);
        memcpy(b, ref->f, sizeof(Tri) * (size_t)ref->nf);
        qsort(a, (size_t)j->ncommitted, sizeof(Tri), tri_cmp);
        qsort(b, (size_t)ref->nf, sizeof(Tri), tri_cmp);
        if (memcmp(a, b, sizeof(Tri) * (size_t)j->ncommitted) != 0)
            return vfail(detail, "weld.face", "committed face set differs from reference");
    }
#endif

    /* normal step: before publish the output buffer must be untouched (all zero) */
    if (!j->normals_published) {
        for (int i = 0; i < j->nu; i++)
            if (j->nrm[i].x != 0.0f || j->nrm[i].y != 0.0f || j->nrm[i].z != 0.0f)
                return vfail(detail, "normal.publish", "nrm[%d] written before publish gate", i);
    } else {
        /* once published, ALL vertices must carry a unit normal or be accounted in bad */
        int zeros = 0;
        for (int i = 0; i < j->nu; i++) {
            float l = len3(j->nrm[i]);
            if (l < 1e-12f) { zeros++; continue; }
            if (fabsf(l - 1.0f) > 1e-5f)
                return vfail(detail, "normal.publish", "nrm[%d] length %.6f neither unit nor zero (half-normalized state)", i, (double)l);
        }
        if (j->stage == ST_DONE) {
            if (zeros != j->bad
#if FAULT_INJECT == 0
                || j->bad != ref->bad
#endif
               )
                return vfail(detail, "normal.publish", "zero-normal count %d inconsistent (bad=%d)", zeros, j->bad);
#if FAULT_INJECT == 0
            for (int i = 0; i < j->nu; i++) {
                V3 d = sub3(j->nrm[i], ref->nrm[i]);
                if (len3(d) > 1e-5f)
                    return vfail(detail, "normal.publish", "nrm[%d] differs from reference", i);
            }
#endif
        } else if (zeros > 0) {
            return vfail(detail, "normal.publish", "%d vertices still zero while gate is on (publish not atomic)", zeros);
        }
    }
    return 1;
}

/* ---------------- test harness ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* 16x16 grid, every vertex duplicated (512 input -> 256 welded), every face
 * duplicated, plus explicit degenerate faces that welding must drop. */
static void build_input(MeshIn* in) {
    enum { N = 16 };
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++) {
            in_add_v(in, (float)i, (float)j, 0.0f);
            in_add_v(in, (float)i, (float)j, 0.0f);   /* exact duplicate */
        }
    for (int j = 0; j < N - 1; j++)
        for (int i = 0; i < N - 1; i++) {
            int p = 2 * (j * N + i);
            int q = p + 2, r = p + 2 * N, s = r + 2;
            in_add_f(in, p, q, s);
            in_add_f(in, p, s, r);
        }
    int base_nf = in->nf;
    for (int k = 0; k < base_nf; k++) in_add_f(in, in->f[k].a, in->f[k].b, in->f[k].c);
    for (int i = 0; i < 15; i++) in_add_f(in, 2 * i, 2 * i, 2 * i + 2);  /* degenerate */
}

#if FAULT_INJECT == 0
static unsigned rng_state;
static unsigned rng_next(void) {
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5;
    return rng_state;
}

typedef struct { int stage, vcur, fcur, acur, ncommitted, nu; } Snapshot;

static void snapshot(const Job* j, Snapshot* s) {
    s->stage = j->stage; s->vcur = j->vcur; s->fcur = j->fcur;
    s->acur = j->acur; s->ncommitted = j->ncommitted; s->nu = j->nu;
}
static int retained(const Job* j, const Snapshot* s) {
    return j->stage >= s->stage && j->vcur >= s->vcur && j->fcur >= s->fcur &&
           j->acur >= s->acur && j->ncommitted >= s->ncommitted && j->nu >= s->nu;
}
#endif

static int run_to_stop(Job* j) {
    for (;;) {
        int r = job_step(j);
        if (r != 0) return r;
    }
}

#if FAULT_INJECT == 0
static int compare_final(const Job* j, const Reference* ref, char* detail) {
    if (j->stage != ST_DONE) { snprintf(detail, 256, "job not done"); return 0; }
    if (!validate(j, ref, detail)) return 0;
    if (!box_equal(&j->box, &ref->box)) { snprintf(detail, 256, "bbox: final box mismatch"); return 0; }
    return 1;
}

/* One trial: cancel at fire1 (and optionally fire2), validating every state the
 * cancellation leaves behind, then resume and require the exact reference. */
static int trial(const MeshIn* in, const Reference* ref, int fire1, int fire2, int use_signal) {
    Job j;
    job_init(&j, in);
    if (use_signal) j.sig_arm = fire1; else j.cancel_request = 1;

    int r = run_to_stop(&j);
    if (r != -1) { printf("  FAIL: expected cancellation at cp %d\n", fire1); return 0; }
    char detail[256];
    if (!validate(&j, ref, detail)) { printf("  FAIL after cancel@%d: %s\n", fire1, detail); return 0; }

    Snapshot before; snapshot(&j, &before);
    if (fire2 >= 0) {
        if (use_signal) j.sig_arm = fire2; else j.cancel_request = 1;
        r = run_to_stop(&j);
        if (r != -1) { printf("  FAIL: expected second cancellation at cp %d\n", fire2); return 0; }
        if (!validate(&j, ref, detail)) { printf("  FAIL after cancel@%d: %s\n", fire2, detail); return 0; }
        Snapshot mid; snapshot(&j, &mid);
        if (!retained(&j, &before)) { printf("  FAIL: second cancel lost committed work\n"); return 0; }
        snapshot(&j, &before); (void)mid;
    }

    /* one step with no cancel pending: committed batches must not be redone */
    r = job_step(&j);
    if (r == -1) { printf("  FAIL: spurious cancellation after resume\n"); return 0; }
    if (r == 0 && !retained(&j, &before)) { printf("  FAIL: resume reverted committed batches\n"); return 0; }

    while (job_step(&j) == 0) { }
    if (!compare_final(&j, ref, detail)) { printf("  FAIL resumed final: %s\n", detail); return 0; }
    return 1;
}
#endif

int main(int argc, char** argv) {
#if FAULT_INJECT != 0
    (void)argc; (void)argv;
#endif
    MeshIn in; memset(&in, 0, sizeof(in));
    build_input(&in);

    Reference ref;
    reference_run(&in, &ref);

    install_handler();
    char detail[256];

    /* regression for the signal() one-shot trap: two deliveries must both be
     * caught; reaching the next line at all proves the handler persisted. */
    raise(SIGINT);
    raise(SIGINT);
    check("handler_persistent", sig_flag == 1);
    sig_flag = 0;

    check("reference_input", in.nv == 512 && in.nf == 915);
    check("reference_weld", ref.nv == 256 && ref.nf == 900 && ref.bad == 0);

#if FAULT_INJECT != 0
    /* Fault demo: sweep every cancellation point and report the first state
     * that violates an invariant, naming the step whose atomicity broke. The
     * clean pipeline crosses well under MAXV boundaries; this fixed bound
     * covers every checkpoint the faulted variants can reach. */
    const char* want = FAULT_INJECT == 1 ? "weld" : FAULT_INJECT == 2 ? "normal" : "bbox";
    int detected = 0;
    for (int cp = 0; cp < MAXV && !detected; cp++) {
        Job j;
        job_init(&j, &in);
        j.fault_armed = 1; j.fault_cp = cp;
        int r = run_to_stop(&j);
        if (r != -1) continue;
        if (!validate(&j, &ref, detail)) {
            if (strncmp(detail, want, strlen(want)) != 0) {
                printf("fault_%d_wrong_step=%s (wanted %s)\n", FAULT_INJECT, detail, want);
                return 1;
            }
            printf("fault_%d_detected_at_cp=%d step=%s\n", FAULT_INJECT, cp, detail);
            detected = 1;
        }
    }
    check("fault_detected", detected);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
#else
    /* uninterrupted run: establishes the answer and the number of batch
     * boundaries a full run crosses. */
    Job base;
    job_init(&base, &in);
    while (job_step(&base) == 0) { }
    int ncp = base.checkpoint_count;
    check("uninterrupted_ok", compare_final(&base, &ref, detail));
    if (failures) { printf("detail: %s\n", detail); return 1; }
    printf("checkpoints_per_run=%d\n", ncp);

    /* 1. deterministic sweep: cancel once at EVERY batch boundary */
    int sweep_ok = 1;
    for (int cp = 0; cp < ncp; cp++) {
        if (!trial(&in, &ref, cp, -1, 0)) { sweep_ok = 0; break; }
    }
    check("sweep_every_boundary", sweep_ok);

    /* 2. 100 random-moment cancellations (including double cancels and real
     * SIGINT delivery), each followed by a full consistency check. */
    rng_state = (unsigned)(argc > 1 ? strtoul(argv[1], NULL, 0) : 0x5eed1234u);
    if (rng_state == 0) rng_state = 1;
    int random_ok = 1;
    int signal_trials = 0;
    for (int t = 0; t < 100; t++) {
        int f1 = (int)(rng_next() % (unsigned)ncp);
        int f2 = -1;
        if ((rng_next() & 1u) && f1 + 2 < ncp) f2 = f1 + 1 + (int)(rng_next() % (unsigned)(ncp - f1 - 1));
        int use_signal = (rng_next() & 1u) != 0;
        if (use_signal) signal_trials++;
        if (!trial(&in, &ref, f1, f2, use_signal)) {
            printf("  in random trial %d (fires %d,%d signal=%d)\n", t, f1, f2, use_signal);
            random_ok = 0;
            break;
        }
    }
    printf("signal_trials=%d\n", signal_trials);
    check("random_100_self_consistent", random_ok);

    /* 3. a cancel flag already pending before the first batch */
    Job j;
    job_init(&j, &in);
    j.cancel_request = 1;
    int r = run_to_stop(&j);
    int early_ok = r == -1 && validate(&j, &ref, detail) && j.checkpoint_count == 1;
    while (job_step(&j) == 0) { }
    early_ok = early_ok && compare_final(&j, &ref, detail);
    check("cancel_before_start", early_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
#endif
}
