/* layout: AoS (struct array) vs SoA (parallel arrays) vertex storage benchmark.
 *
 * Vertex = 32 B: position(3f) + normal(3f) + uv(2f).
 *   AoS: Vertex a[n]
 *   SoA: V3 pos[n]; V3 nrm[n]; V2 uv[n];   (three parallel arrays, GPU-stream style)
 *
 * Access patterns, each on both layouts:
 *   trav_full   sequential, all fields, sum        (transform-like)
 *   trav_pos    sequential, position only, bbox    (compute-heavy traversal)
 *   trav_x      sequential, pos.x only, sum        (bandwidth probe)
 *   rand_full   random order (independent), one field per group
 *   rand_pos    random order (independent), position only
 *   chase_pos   random order (dependent chain), pos.x -> next index
 *   chase_full  dependent chain, next = f(pos.x, nrm.x, uv.u): all 3 lines needed
 *
 * Reported: ns per vertex-visit (trav/rand) or ns per step (chase), min of batches.
 * Sizes above 100k are out-of-spec reference points to expose the DRAM asymptote.
 */
#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef struct { float x, y, z; } V3;
typedef struct { float u, v; } V2;
typedef struct { V3 pos; V3 nrm; V2 uv; } Vertex;   /* 32 B, two per 64 B line */

typedef struct { V3* pos; V3* nrm; V2* uv; } SoA;

static volatile float g_sink;

static uint32_t xs = 0x9e3779b9u;
static uint32_t xr(void) { xs ^= xs << 13; xs ^= xs >> 17; xs ^= xs << 5; return xs; }
static float frand(void) { return (float)(xr() >> 8) * (1.0f / 16777216.0f); }

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void* xalloc(size_t bytes) {
    void* p = NULL;
    if (posix_memalign(&p, 64, (bytes + 63) & ~(size_t)63) != 0) { fprintf(stderr, "oom\n"); exit(1); }
    return p;
}

static uint32_t f2u(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }
static float u2f(uint32_t u) { float f; memcpy(&f, &u, 4); return f; }

/* ---------------- traversal kernels ---------------- */

static float trav_full_aos(const void* d, size_t n) {
    const Vertex* a = d;
    float s0 = 0, s1 = 0, s2 = 0;
    for (size_t i = 0; i < n; i++) {
        s0 += a[i].pos.x + a[i].pos.y + a[i].pos.z;
        s1 += a[i].nrm.x + a[i].nrm.y + a[i].nrm.z;
        s2 += a[i].uv.u + a[i].uv.v;
    }
    return s0 + s1 + s2;
}

static float trav_full_soa(const void* d, size_t n) {
    const SoA* s = d;
    float s0 = 0, s1 = 0, s2 = 0;
    for (size_t i = 0; i < n; i++) {
        s0 += s->pos[i].x + s->pos[i].y + s->pos[i].z;
        s1 += s->nrm[i].x + s->nrm[i].y + s->nrm[i].z;
        s2 += s->uv[i].u + s->uv[i].v;
    }
    return s0 + s1 + s2;
}

static float trav_pos_aos(const void* d, size_t n) {   /* bbox grow */
    const Vertex* a = d;
    float mnx = 1e30f, mny = 1e30f, mnz = 1e30f;
    float mxx = -1e30f, mxy = -1e30f, mxz = -1e30f;
    for (size_t i = 0; i < n; i++) {
        float x = a[i].pos.x, y = a[i].pos.y, z = a[i].pos.z;
        mnx = mnx < x ? mnx : x;  mxx = mxx > x ? mxx : x;
        mny = mny < y ? mny : y;  mxy = mxy > y ? mxy : y;
        mnz = mnz < z ? mnz : z;  mxz = mxz > z ? mxz : z;
    }
    return mnx + mny + mnz + mxx + mxy + mxz;
}

static float trav_pos_soa(const void* d, size_t n) {
    const SoA* s = d;
    float mnx = 1e30f, mny = 1e30f, mnz = 1e30f;
    float mxx = -1e30f, mxy = -1e30f, mxz = -1e30f;
    for (size_t i = 0; i < n; i++) {
        float x = s->pos[i].x, y = s->pos[i].y, z = s->pos[i].z;
        mnx = mnx < x ? mnx : x;  mxx = mxx > x ? mxx : x;
        mny = mny < y ? mny : y;  mxy = mxy > y ? mxy : y;
        mnz = mnz < z ? mnz : z;  mxz = mxz > z ? mxz : z;
    }
    return mnx + mny + mnz + mxx + mxy + mxz;
}

static float trav_x_aos(const void* d, size_t n) {   /* lightest possible touch */
    const Vertex* a = d;
    float s = 0;
    for (size_t i = 0; i < n; i++) s += a[i].pos.x;
    return s;
}

static float trav_x_soa(const void* d, size_t n) {
    const SoA* s = d;
    float r = 0;
    for (size_t i = 0; i < n; i++) r += s->pos[i].x;
    return r;
}

/* ---------------- random, independent (throughput) ---------------- */

static float rand_full_aos(const void* d, const uint32_t* p, size_t n) {
    const Vertex* a = d;
    float s = 0;
    for (size_t i = 0; i < n; i++) {
        const Vertex* v = &a[p[i]];
        s += v->pos.x + v->nrm.x + v->uv.u;   /* one field per group: 1 cache line */
    }
    return s;
}

static float rand_full_soa(const void* d, const uint32_t* p, size_t n) {
    const SoA* s = d;
    float r = 0;
    for (size_t i = 0; i < n; i++) {
        uint32_t j = p[i];
        r += s->pos[j].x + s->nrm[j].x + s->uv[j].u;   /* 3 separate cache lines */
    }
    return r;
}

static float rand_pos_aos(const void* d, const uint32_t* p, size_t n) {
    const Vertex* a = d;
    float s = 0;
    for (size_t i = 0; i < n; i++) {
        const V3* q = &a[p[i]].pos;
        s += q->x + q->y + q->z;
    }
    return s;
}

static float rand_pos_soa(const void* d, const uint32_t* p, size_t n) {
    const SoA* s = d;
    float r = 0;
    for (size_t i = 0; i < n; i++) {
        const V3* q = &s->pos[p[i]];
        r += q->x + q->y + q->z;
    }
    return r;
}

/* ---------------- random, dependent chain (latency) ----------------
 * next index is bit-encoded into the vertex fields, so each step's loads
 * cannot start until the previous step's loads complete: the prefetcher
 * and OOO window cannot hide anything. Pure cache-level latency.
 */
static float chase_pos_aos(const void* d, size_t steps) {
    const Vertex* a = d;
    uint32_t cur = 0;
    float acc = 0;
    for (size_t i = 0; i < steps; i++) {
        const Vertex* v = &a[cur];
        cur = f2u(v->pos.x);
        acc += v->pos.y;
    }
    return acc;
}

static float chase_pos_soa(const void* d, size_t steps) {
    const SoA* s = d;
    uint32_t cur = 0;
    float acc = 0;
    for (size_t i = 0; i < steps; i++) {
        cur = f2u(s->pos[cur].x);
        acc += s->pos[cur].y;   /* note: reads next element's y, still 1 line/step */
    }
    return acc;
}

static float chase_full_aos(const void* d, size_t steps) {
    const Vertex* a = d;
    uint32_t cur = 0;
    float acc = 0;
    for (size_t i = 0; i < steps; i++) {
        const Vertex* v = &a[cur];
        cur = f2u(v->pos.x) ^ f2u(v->nrm.x) ^ f2u(v->uv.u);   /* all in the same 32 B */
        acc += v->pos.y;
    }
    return acc;
}

static float chase_full_soa(const void* d, size_t steps) {
    const SoA* s = d;
    uint32_t cur = 0;
    float acc = 0;
    for (size_t i = 0; i < steps; i++) {
        uint32_t j = cur;
        cur = f2u(s->pos[j].x) ^ f2u(s->nrm[j].x) ^ f2u(s->uv[j].u);   /* 3 lines */
        acc += s->pos[j].y;
    }
    return acc;
}

/* ---------------- harness ---------------- */

typedef float (*TravKernel)(const void*, size_t);
typedef float (*RandKernel)(const void*, const uint32_t*, size_t);

static double bench_trav(TravKernel k, const void* d, size_t n) {
    size_t reps = 20000000 / n; if (reps < 1) reps = 1;
    double best = 1e30;
    for (int b = 0; b < 7; b++) {
        double t0 = now_s();
        float acc = 0;
        for (size_t r = 0; r < reps; r++) acc += k(d, n);
        double dt = (now_s() - t0) / (double)reps;
        g_sink += acc;
        if (dt < best) best = dt;
    }
    return best / (double)n * 1e9;   /* ns per vertex */
}

static double bench_rand(RandKernel k, const void* d, const uint32_t* p, size_t n) {
    size_t reps = 2000000 / n; if (reps < 1) reps = 1;
    double best = 1e30;
    for (int b = 0; b < 7; b++) {
        double t0 = now_s();
        float acc = 0;
        for (size_t r = 0; r < reps; r++) acc += k(d, p, n);
        double dt = (now_s() - t0) / (double)reps;
        g_sink += acc;
        if (dt < best) best = dt;
    }
    return best / (double)n * 1e9;
}

typedef float (*ChaseKernel)(const void*, size_t);

static double bench_chase(ChaseKernel k, const void* d, size_t steps) {
    double best = 1e30;
    for (int b = 0; b < 5; b++) {
        double t0 = now_s();
        g_sink += k(d, steps);
        double dt = now_s() - t0;
        if (dt < best) best = dt;
    }
    return best / (double)steps * 1e9;   /* ns per step */
}

int main(void) {
    static const size_t sizes[] = { 100, 300, 1000, 3000, 10000, 30000, 100000,
                                    300000, 1000000, 3000000 };   /* last 3: beyond spec */
    const size_t CHASE_STEPS = 1000000;
    printf("sizeof(Vertex)=%zu  L1d=48K(1536v) L2=1M(32kv) L3=32M(1Mv)\n", sizeof(Vertex));

    printf("\n[traversal, ns/vertex]  full=all fields, pos=bbox, x=pos.x only\n");
    printf("%9s | %5s %5s | %5s %5s | %5s %5s\n", "n", "fA", "fS", "pA", "pS", "xA", "xS");
    printf("----------+-------------+-------------+-------------\n");

    /* tables are printed per-section; we compute all sections per size in one pass */
    struct row { size_t n; double v[16]; } rows[16];
    size_t nr = 0;

    for (size_t si = 0; si < sizeof(sizes) / sizeof(sizes[0]); si++) {
        size_t n = sizes[si];
        Vertex* a = xalloc(n * sizeof(Vertex));
        SoA s;
        s.pos = xalloc(n * sizeof(V3));
        s.nrm = xalloc(n * sizeof(V3));
        s.uv  = xalloc(n * sizeof(V2));
        uint32_t* perm = xalloc(n * sizeof(uint32_t));

        for (size_t i = 0; i < n; i++) {
            V3 p = { frand(), frand(), frand() };
            V3 m = { frand(), frand(), frand() };
            V2 t = { frand(), frand() };
            a[i].pos = p; a[i].nrm = m; a[i].uv = t;
            s.pos[i] = p; s.nrm[i] = m; s.uv[i] = t;
        }
        for (size_t i = 0; i < n; i++) perm[i] = (uint32_t)i;
        for (size_t i = n - 1; i > 0; i--) {
            size_t j = xr() % (i + 1);
            uint32_t tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
        }

        struct row* R = &rows[nr++];
        R->n = n;
        R->v[0] = bench_trav(trav_full_aos, a, n);
        R->v[1] = bench_trav(trav_full_soa, &s, n);
        R->v[2] = bench_trav(trav_pos_aos, a, n);
        R->v[3] = bench_trav(trav_pos_soa, &s, n);
        R->v[4] = bench_trav(trav_x_aos, a, n);
        R->v[5] = bench_trav(trav_x_soa, &s, n);
        R->v[6] = bench_rand(rand_full_aos, a, perm, n);
        R->v[7] = bench_rand(rand_full_soa, &s, perm, n);
        R->v[8] = bench_rand(rand_pos_aos, a, perm, n);
        R->v[9] = bench_rand(rand_pos_soa, &s, perm, n);

        /* encode chase chain: pos.x = next index; nrm.x/uv.u are read by
         * chase_full too (forcing their cache lines in), but decode to 0 */
        for (size_t i = 0; i < n; i++) {
            a[i].pos.x = u2f(perm[i]); a[i].nrm.x = u2f(0); a[i].uv.u = u2f(0);
            s.pos[i].x = u2f(perm[i]); s.nrm[i].x = u2f(0); s.uv[i].u = u2f(0);
        }
        R->v[10] = bench_chase((ChaseKernel)chase_pos_aos, a, CHASE_STEPS);
        R->v[11] = bench_chase((ChaseKernel)chase_pos_soa, &s, CHASE_STEPS);
        R->v[12] = bench_chase((ChaseKernel)chase_full_aos, a, CHASE_STEPS);
        R->v[13] = bench_chase((ChaseKernel)chase_full_soa, &s, CHASE_STEPS);

        free(a); free(s.pos); free(s.nrm); free(s.uv); free(perm);
    }

    for (size_t i = 0; i < nr; i++)
        printf("%9zu | %5.2f %5.2f | %5.2f %5.2f | %5.2f %5.2f\n",
               rows[i].n, rows[i].v[0], rows[i].v[1], rows[i].v[2], rows[i].v[3], rows[i].v[4], rows[i].v[5]);

    printf("\n[random independent, ns/visit]  full=1 field per group, pos=position only\n");
    printf("%9s | %5s %5s | %5s %5s\n", "n", "fA", "fS", "pA", "pS");
    printf("----------+-------------+-------------\n");
    for (size_t i = 0; i < nr; i++)
        printf("%9zu | %5.2f %5.2f | %5.2f %5.2f\n",
               rows[i].n, rows[i].v[6], rows[i].v[7], rows[i].v[8], rows[i].v[9]);

    printf("\n[random dependent chain, ns/step]  pos=chain via pos.x, full=chain needs all 3 groups\n");
    printf("%9s | %6s %6s | %6s %6s\n", "n", "pA", "pS", "fA", "fS");
    printf("----------+---------------+---------------\n");
    for (size_t i = 0; i < nr; i++)
        printf("%9zu | %6.2f %6.2f | %6.2f %6.2f\n",
               rows[i].n, rows[i].v[10], rows[i].v[11], rows[i].v[12], rows[i].v[13]);

    return 0;
}
