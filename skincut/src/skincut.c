
/* skincut: influence-ranked bone pruning for skeleton export.
 *
 * influence(b) = weight_sum(b) * (weight_sum(b) / vertex_count(b))
 *              = total weight x average weight per influenced vertex,
 * i.e. decided jointly by influenced-vertex count and weight sum.
 * Prune lowest influence first, remap cut bones' skin weights to the
 * nearest surviving ancestor, then measure animation error vs the full
 * skeleton. If finger bones (few verts, full weight) get pruned first,
 * the ranking used vertex count alone -- the report lists verts and
 * weight_sum side by side so the cause is visible.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAXB 256
#define MAXV 4608
#define MAXW 4
#define FRAMES 16
#define EXPORT_LIMIT 43

typedef struct { float x, y, z; } V3;
typedef struct { float m[4][4]; } M4;

typedef struct {
    char name[32];
    int parent;             /* -1 = root */
    V3 rest_t;              /* local rest translation */
    V3 axis;                /* animation rotation axis */
    float amp, freq, phase; /* angle = amp * sin(2*pi*freq*t + phase) */
    int protect;            /* never pruned */
} Bone;

typedef struct { int bone[MAXW]; float w[MAXW]; int n; } Skin;

typedef struct {
    char name[32];
    Bone b[MAXB]; int nb;
    V3 v[MAXV]; Skin s[MAXV]; int nv;
} Model;

/* ---- 4x4 rigid transforms ---- */

static M4 m4_ident(void) {
    M4 r; memset(&r, 0, sizeof r);
    for (int i = 0; i < 4; i++) r.m[i][i] = 1.0f;
    return r;
}

static M4 m4_trans(V3 t) {
    M4 r = m4_ident();
    r.m[0][3] = t.x; r.m[1][3] = t.y; r.m[2][3] = t.z;
    return r;
}

static M4 m4_rot(V3 a, float ang) {
    float l = sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);
    if (l < 1e-9f) return m4_ident();
    float x = a.x/l, y = a.y/l, z = a.z/l;
    float c = cosf(ang), s = sinf(ang), t = 1.0f - c;
    M4 r = m4_ident();
    r.m[0][0] = t*x*x + c;   r.m[0][1] = t*x*y - s*z; r.m[0][2] = t*x*z + s*y;
    r.m[1][0] = t*x*y + s*z; r.m[1][1] = t*y*y + c;   r.m[1][2] = t*y*z - s*x;
    r.m[2][0] = t*x*z - s*y; r.m[2][1] = t*y*z + s*x; r.m[2][2] = t*z*z + c;
    return r;
}

static M4 m4_mul(M4 a, M4 b) {
    M4 r;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            r.m[i][j] = 0.0f;
            for (int k = 0; k < 4; k++) r.m[i][j] += a.m[i][k] * b.m[k][j];
        }
    return r;
}

static V3 m4_apply(M4 m, V3 p) {
    V3 r;
    r.x = m.m[0][0]*p.x + m.m[0][1]*p.y + m.m[0][2]*p.z + m.m[0][3];
    r.y = m.m[1][0]*p.x + m.m[1][1]*p.y + m.m[1][2]*p.z + m.m[1][3];
    r.z = m.m[2][0]*p.x + m.m[2][1]*p.y + m.m[2][2]*p.z + m.m[2][3];
    return r;
}

static M4 m4_inv_rigid(M4 m) {
    M4 r = m4_ident();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            r.m[i][j] = m.m[j][i];
    for (int i = 0; i < 3; i++)
        r.m[i][3] = -(r.m[i][0]*m.m[0][3] + r.m[i][1]*m.m[1][3] + r.m[i][2]*m.m[2][3]);
    return r;
}

/* ---- model construction (deterministic, no rand) ---- */

static int add_bone(Model* m, const char* name, int parent, V3 t,
                    V3 axis, float amp, float freq, float phase, int protect) {
    if (m->nb >= MAXB) return -1;
    Bone* b = &m->b[m->nb];
    snprintf(b->name, sizeof b->name, "%s", name);
    b->parent = parent; b->rest_t = t;
    b->axis = axis; b->amp = amp; b->freq = freq; b->phase = phase;
    b->protect = protect;
    return m->nb++;
}

static int add_vert(Model* m, V3 p, const int* bones, const float* ws, int n) {
    if (m->nv >= MAXV || n > MAXW) return -1;
    m->v[m->nv] = p;
    Skin* s = &m->s[m->nv];
    s->n = n;
    for (int i = 0; i < n; i++) { s->bone[i] = bones[i]; s->w[i] = ws[i]; }
    return m->nv++;
}

static void pose_rest(const Model* m, M4* g) {
    for (int b = 0; b < m->nb; b++) {
        M4 local = m4_trans(m->b[b].rest_t);
        g[b] = m->b[b].parent < 0 ? local : m4_mul(g[m->b[b].parent], local);
    }
}

static void pose_anim(const Model* m, float t, M4* g) {
    for (int b = 0; b < m->nb; b++) {
        const Bone* bn = &m->b[b];
        float ang = bn->amp * sinf(2.0f * (float)M_PI * bn->freq * t + bn->phase);
        M4 local = m4_mul(m4_trans(bn->rest_t), m4_rot(bn->axis, ang));
        g[b] = bn->parent < 0 ? local : m4_mul(g[bn->parent], local);
    }
}

static V3 jitter(int bone, int k) {
    V3 r = { 0.03f*sinf(bone*1.7f + k*2.1f),
             0.03f*cosf(bone*0.9f + k*1.3f),
             0.03f*sinf(bone*0.5f + k*0.7f) };
    return r;
}

/* flesh out one bone with vertices weighted to it (and optionally its parent) */
static void flesh(Model* m, const M4* g_rest, int bi, int nverts, float self_w) {
    V3 base = { g_rest[bi].m[0][3], g_rest[bi].m[1][3], g_rest[bi].m[2][3] };
    int parent = m->b[bi].parent;
    for (int k = 0; k < nverts; k++) {
        V3 p = base;
        V3 j = jitter(bi, k);
        p.x += j.x; p.y += j.y; p.z += j.z;
        if (self_w >= 1.0f || parent < 0) {
            int wb[1] = { bi }; float ww[1] = { 1.0f };
            add_vert(m, p, wb, ww, 1);
        } else {
            int wb[2] = { bi, parent }; float ww[2] = { self_w, 1.0f - self_w };
            add_vert(m, p, wb, ww, 2);
        }
    }
}

/* humanoid: root + 4 spine + 2x(shoulder,upper,lower,hand + 5 fingers x 3 joints)
 * = 43 structural bones, plus nhelper twist/cloth helper bones.
 * Fingers: few vertices, weight 1.0. Helpers: many vertices, weight 0.01. */
static void build_model(Model* m, const char* name, int nhelper) {
    memset(m, 0, sizeof *m);
    snprintf(m->name, sizeof m->name, "%s", name);
    V3 Z = {0,0,1}, Y = {0,1,0}, HX = {0.3f,0.5f,0.8f};

    add_bone(m, "root", -1, (V3){0,0,0}, Z, 0.0f, 0.0f, 0.0f, 1);
    for (int i = 0; i < 4; i++) {
        char nm[32]; snprintf(nm, sizeof nm, "spine_%d", i);
        add_bone(m, nm, i == 0 ? 0 : i, (V3){0, 0.25f, 0}, Z, 0.15f, 1.0f, i*0.4f, 0);
    }
    int sh[2];
    for (int side = 0; side < 2; side++) {
        float sx = side == 0 ? 1.0f : -1.0f;
        const char* tag = side == 0 ? "L" : "R";
        char nm[32];
        snprintf(nm, sizeof nm, "shoulder_%s", tag);
        sh[side] = add_bone(m, nm, 4, (V3){sx*0.18f, 0.02f, 0}, Z, 0.25f, 1.0f, 0.3f, 0);
        snprintf(nm, sizeof nm, "upperarm_%s", tag);
        int ua = add_bone(m, nm, sh[side], (V3){sx*0.22f, 0, 0}, Z, 0.25f, 1.0f, 0.6f, 0);
        snprintf(nm, sizeof nm, "lowerarm_%s", tag);
        int la = add_bone(m, nm, ua, (V3){sx*0.20f, 0, 0}, Z, 0.25f, 1.0f, 0.9f, 0);
        snprintf(nm, sizeof nm, "hand_%s", tag);
        int hand = add_bone(m, nm, la, (V3){sx*0.16f, 0, 0}, Z, 0.25f, 1.0f, 1.2f, 0);
        int prev[5];
        for (int f = 0; f < 5; f++) {
            prev[f] = hand;
            for (int j = 0; j < 3; j++) {
                snprintf(nm, sizeof nm, "finger_%s_%d_%d", tag, f, j);
                V3 t = j == 0 ? (V3){sx*0.05f, (f-2)*0.035f, 0.01f*(f%2)}
                              : (V3){sx*0.035f, 0, 0};
                prev[f] = add_bone(m, nm, prev[f], t, Y, 0.6f, 2.0f, f*0.5f + j*0.25f, 0);
            }
        }
    }
    int hpar[6] = { 1, 2, 3, 4, sh[0], sh[1] };
    for (int i = 0; i < nhelper; i++) {
        char nm[32]; snprintf(nm, sizeof nm, "helper_%d", i);
        V3 t = { 0.05f*sinf(i*1.3f), 0.04f*cosf(i*0.7f), 0.05f*sinf(i*0.9f) };
        add_bone(m, nm, hpar[i % 6], t, HX, 0.8f, 1.5f, i*0.7f, 0);
    }

    M4 g_rest[MAXB];
    pose_rest(m, g_rest);
    for (int b = 1; b <= 4; b++) flesh(m, g_rest, b, 24, 0.5f);        /* spine: blended */
    for (int b = 5; b < 43; b++) {
        if (strstr(m->b[b].name, "finger")) flesh(m, g_rest, b, 6, 1.0f);  /* rigid */
        else                                flesh(m, g_rest, b, 16, 0.5f); /* arm chain */
    }
    for (int b = 43; b < m->nb; b++) flesh(m, g_rest, b, 24, 0.01f);   /* helpers: tiny weight */
}

/* ---- influence ---- */

static void calc_influence(const Model* m, int* cnt, double* ws, double* infl) {
    for (int b = 0; b < m->nb; b++) { cnt[b] = 0; ws[b] = 0.0; }
    for (int v = 0; v < m->nv; v++)
        for (int i = 0; i < m->s[v].n; i++) {
            cnt[m->s[v].bone[i]]++;
            ws[m->s[v].bone[i]] += m->s[v].w[i];
        }
    for (int b = 0; b < m->nb; b++)
        infl[b] = cnt[b] > 0 ? ws[b] * (ws[b] / cnt[b]) : 0.0;
}

static double g_key[MAXB];
static int cmp_idx(const void* a, const void* b) {
    double d = g_key[*(const int*)a] - g_key[*(const int*)b];
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

/* pick the nb-limit lowest-influence bones; cutlist ordered by influence asc.
 * count_only=1 reproduces the buggy metric (vertex count alone). */
static int prune_bones(const Model* m, int limit, int count_only,
                       const int* cnt, const double* infl, int* cutlist) {
    int idx[MAXB], n = 0;
    for (int b = 0; b < m->nb; b++) {
        if (m->b[b].protect) continue;
        g_key[b] = count_only ? (double)cnt[b] : infl[b];
        idx[n++] = b;
    }
    qsort(idx, n, sizeof(int), cmp_idx);
    int ncut = m->nb - limit;
    if (ncut < 0) ncut = 0;
    if (ncut > n) ncut = n;
    for (int i = 0; i < ncut; i++) cutlist[i] = idx[i];
    return ncut;
}

static void build_keep_map(const Model* m, const int* cutlist, int ncut, int* keep_map) {
    int cut[MAXB] = {0};
    for (int i = 0; i < ncut; i++) cut[cutlist[i]] = 1;
    for (int b = 0; b < m->nb; b++) {
        int k = b;
        while (k >= 0 && cut[k]) k = m->b[k].parent;
        keep_map[b] = k < 0 ? 0 : k;
    }
}

/* remap skin weights through keep_map and renormalize */
static void remap_skins(const Model* m, const int* keep_map, Skin* out) {
    for (int v = 0; v < m->nv; v++) {
        Skin* d = &out[v];
        d->n = 0;
        float sum = 0.0f;
        for (int i = 0; i < m->s[v].n; i++) {
            int kb = keep_map[m->s[v].bone[i]];
            float w = m->s[v].w[i];
            int j;
            for (j = 0; j < d->n; j++) if (d->bone[j] == kb) break;
            if (j < d->n) d->w[j] += w;
            else { d->bone[j] = kb; d->w[j] = w; d->n++; }
            sum += w;
        }
        if (sum > 1e-9f)
            for (int i = 0; i < d->n; i++) d->w[i] /= sum;
    }
}

static V3 skin_vert(const M4* S, V3 p, const Skin* s) {
    V3 r = {0,0,0};
    for (int i = 0; i < s->n; i++) {
        V3 q = m4_apply(S[s->bone[i]], p);
        r.x += s->w[i]*q.x; r.y += s->w[i]*q.y; r.z += s->w[i]*q.z;
    }
    return r;
}

/* RMS vertex displacement between original and remapped skins over the
 * animation, normalized by the rest-pose bbox diagonal. */
static double anim_error(const Model* m, const Skin* alt, double* max_out) {
    M4 g_rest[MAXB], inv_rest[MAXB], g[MAXB], S[MAXB];
    pose_rest(m, g_rest);
    for (int b = 0; b < m->nb; b++) inv_rest[b] = m4_inv_rigid(g_rest[b]);

    V3 lo = m->v[0], hi = m->v[0];
    for (int v = 0; v < m->nv; v++) {
        if (m->v[v].x < lo.x) lo.x = m->v[v].x;
        if (m->v[v].x > hi.x) hi.x = m->v[v].x;
        if (m->v[v].y < lo.y) lo.y = m->v[v].y;
        if (m->v[v].y > hi.y) hi.y = m->v[v].y;
        if (m->v[v].z < lo.z) lo.z = m->v[v].z;
        if (m->v[v].z > hi.z) hi.z = m->v[v].z;
    }
    double diag = sqrt((hi.x-lo.x)*(hi.x-lo.x) + (hi.y-lo.y)*(hi.y-lo.y) + (hi.z-lo.z)*(hi.z-lo.z));
    if (diag < 1e-9) diag = 1.0;

    double acc = 0.0, maxd = 0.0;
    long n = 0;
    for (int f = 0; f < FRAMES; f++) {
        pose_anim(m, f / (float)FRAMES, g);
        for (int b = 0; b < m->nb; b++) S[b] = m4_mul(g[b], inv_rest[b]);
        for (int v = 0; v < m->nv; v++) {
            V3 p0 = skin_vert(S, m->v[v], &m->s[v]);
            V3 p1 = skin_vert(S, m->v[v], &alt[v]);
            double dx = p0.x-p1.x, dy = p0.y-p1.y, dz = p0.z-p1.z;
            double d2 = dx*dx + dy*dy + dz*dz;
            acc += d2; n++;
            if (d2 > maxd) maxd = d2;
        }
    }
    if (max_out) *max_out = sqrt(maxd) / diag;
    return sqrt(acc / n) / diag;
}

/* list cut bones with verts and weight_sum side by side; returns 1 if any
 * finger bone was pruned (the signature of a vertex-count-only ranking). */
static int report_cuts(const Model* m, const int* cutlist, int ncut,
                       const int* cnt, const double* ws, const double* infl) {
    int finger_cut = 0;
    for (int i = 0; i < ncut; i++) {
        int b = cutlist[i];
        printf("cut_bone name=%-14s verts=%3d weight_sum=%8.4f influence=%9.6f\n",
               m->b[b].name, cnt[b], ws[b], infl[b]);
        if (strstr(m->b[b].name, "finger")) finger_cut = 1;
    }
    if (finger_cut)
        printf("WARNING=finger bones pruned -- influence likely ranked by vertex "
               "count only; compare verts vs weight_sum above\n");
    return finger_cut;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static Model g_model;
static Skin g_remapped[MAXV];

int main(void) {
    const char* names[3] = { "imp", "soldier", "dragon" };
    int helpers[3] = { 8, 60, 160 };            /* bone counts 51 / 103 / 203 */
    double result[3][3];                        /* cut_count, cut_influence_sum, anim_error */

    for (int mi = 0; mi < 3; mi++) {
        Model* m = &g_model;
        build_model(m, names[mi], helpers[mi]);

        int cnt[MAXB]; double ws[MAXB], infl[MAXB];
        calc_influence(m, cnt, ws, infl);

        int cutlist[MAXB];
        int ncut = prune_bones(m, EXPORT_LIMIT, 0, cnt, infl, cutlist);
        int keep_map[MAXB];
        build_keep_map(m, cutlist, ncut, keep_map);
        remap_skins(m, keep_map, g_remapped);
        double maxe, err = anim_error(m, g_remapped, &maxe);

        double infl_sum = 0.0;
        for (int i = 0; i < ncut; i++) infl_sum += infl[cutlist[i]];

        printf("model=%s bones=%d limit=%d\n", m->name, m->nb, EXPORT_LIMIT);
        printf("cut_count=%d\n", ncut);
        printf("cut_influence_sum=%.6f\n", infl_sum);
        printf("anim_error=%.6f\n", err);
        int finger_cut = report_cuts(m, cutlist, ncut, cnt, ws, infl);

        char ck[64];
        snprintf(ck, sizeof ck, "%s_cut_count", m->name);
        check(ck, ncut == helpers[mi]);
        snprintf(ck, sizeof ck, "%s_no_finger_cut", m->name);
        check(ck, !finger_cut);
        snprintf(ck, sizeof ck, "%s_error_small", m->name);
        check(ck, err < 0.02);

        result[mi][0] = ncut; result[mi][1] = infl_sum; result[mi][2] = err;
    }

    /* control: the buggy metric (vertex count only) must trip the diagnostic */
    {
        Model* m = &g_model;
        build_model(m, "soldier", 60);
        int cnt[MAXB]; double ws[MAXB], infl[MAXB];
        calc_influence(m, cnt, ws, infl);
        int cutlist[MAXB];
        int ncut = prune_bones(m, EXPORT_LIMIT, 1, cnt, infl, cutlist);
        int keep_map[MAXB];
        build_keep_map(m, cutlist, ncut, keep_map);
        remap_skins(m, keep_map, g_remapped);
        double err = anim_error(m, g_remapped, 0);
        printf("control=count_only_metric\n");
        int finger_cut = report_cuts(m, cutlist, ncut, cnt, ws, infl);
        check("countonly_trips_finger_warning", finger_cut);
        check("countonly_error_much_worse", err > 10.0 * result[1][2]);
    }

    printf("== summary: model cut_count cut_influence_sum anim_error ==\n");
    for (int mi = 0; mi < 3; mi++)
        printf("%s %d %.6f %.6f\n", names[mi],
               (int)result[mi][0], result[mi][1], result[mi][2]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
