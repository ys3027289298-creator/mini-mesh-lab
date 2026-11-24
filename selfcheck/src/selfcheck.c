/* selfcheck: post-import degeneracy screen.
 *
 * Four classes of import corruption, each with a derived threshold:
 *
 *  1. zero vertices     finite vertex count < 3.
 *                       3 = minimum to form even one face; not a tuned knob.
 *  2. zero bbox volume  min extent / diagonal < 1e-4.
 *                       float32 import-chain noise floors out near 1e-6
 *                       (2^-23 relative, a few multiply-adds); the thinnest
 *                       legitimate parts (sheet metal ~0.1mm on 100mm) sit
 *                       near 1e-3. 1e-4 splits the three-decade gap.
 *  3. normals aligned   mean resultant length R = |sum(unit n)| / N > 0.95.
 *                       R = 1 means all identical; R = 0.95 corresponds to
 *                       normals confined to a ~26 deg cone (R = (1+cos a)/2),
 *                       i.e. carrying no shape information for a solid.
 *  4. isolated vertices unreferenced vertex farther than s/2 from every face.
 *                       s = sqrt(2A / (sqrt(3) V)) comes from vertex density:
 *                       closed mesh Euler gives F ~= 2V, mean triangle area
 *                       A/2V, equilateral edge for that area is s. A stray
 *                       vertex is either on the surface (d ~ sagitta ~ 0) or
 *                       off by a lattice step (d ~ s); s/2 is the midpoint
 *                       between the two shells, maximal margin on both sides.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 400
#define MAXF 600
#define PI 3.14159265358979f

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct {
    V3 v[MAXV]; int nv;
    Tri f[MAXF]; int nf;
    V3 n[MAXV]; int has_normals;
} Mesh;

static V3 sub3(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 cross3(V3 p, V3 q) { V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r; }
static float dot3(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float len3(V3 p) { return sqrtf(dot3(p, p)); }
static int finite3(V3 p) { return isfinite(p.x) && isfinite(p.y) && isfinite(p.z); }

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

static int add_v(Mesh* m, float x, float y, float z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv].x = x; m->v[m->nv].y = y; m->v[m->nv].z = z;
    return m->nv++;
}

static int add_f(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a<0||b<0||c<0||a>=m->nv||b>=m->nv||c>=m->nv) return -1;
    if (a==b||b==c||a==c) return -1;
    m->f[m->nf].a=a; m->f[m->nf].b=b; m->f[m->nf].c=c;
    return m->nf++;
}

/* ---------------- class metrics (threshold-free) ---------------- */

/* class 1: vertices that survived parsing as finite floats */
static int finite_verts(const Mesh* m) {
    int c = 0;
    for (int i = 0; i < m->nv; i++)
        if (finite3(m->v[i])) c++;
    return c;
}

/* class 2: min bbox extent / diagonal; 0 when there is no usable geometry */
static float bbox_ratio(const Mesh* m) {
    V3 lo = { 1e30f, 1e30f, 1e30f }, hi = { -1e30f, -1e30f, -1e30f };
    int n = 0;
    for (int i = 0; i < m->nv; i++) {
        V3 p = m->v[i];
        if (!finite3(p)) continue;
        if (p.x < lo.x) lo.x = p.x;
        if (p.x > hi.x) hi.x = p.x;
        if (p.y < lo.y) lo.y = p.y;
        if (p.y > hi.y) hi.y = p.y;
        if (p.z < lo.z) lo.z = p.z;
        if (p.z > hi.z) hi.z = p.z;
        n++;
    }
    if (!n) return 0.0f;
    V3 s = sub3(hi, lo);
    float d = len3(s);
    if (d <= 0.0f) return 0.0f;
    float mn = s.x < s.y ? s.x : s.y;
    if (s.z < mn) mn = s.z;
    return mn / d;
}

/* class 3: mean resultant length of vertex normals; -1 when not applicable */
static float normal_r(const Mesh* m) {
    if (!m->has_normals) return -1.0f;
    double sx = 0, sy = 0, sz = 0;
    int c = 0;
    for (int i = 0; i < m->nv; i++) {
        V3 n = m->n[i];
        if (!finite3(n)) continue;
        float l = len3(n);
        if (l < 1e-12f) continue;
        sx += n.x / l; sy += n.y / l; sz += n.z / l;
        c++;
    }
    if (!c) return -1.0f;
    return (float)(sqrt(sx*sx + sy*sy + sz*sz) / c);
}

/* squared distance point to segment */
static float d2_seg(V3 p, V3 a, V3 b) {
    V3 ab = sub3(b, a), ap = sub3(p, a);
    float d = dot3(ab, ab);
    float t = d > 0.0f ? dot3(ap, ab) / d : 0.0f;
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    V3 q = { a.x + t*ab.x, a.y + t*ab.y, a.z + t*ab.z };
    return dot3(sub3(p, q), sub3(p, q));
}

/* squared distance point to triangle (Ericson, Real-Time Collision Detection) */
static float d2_tri(V3 p, V3 a, V3 b, V3 c) {
    V3 ab = sub3(b, a), ac = sub3(c, a), ap = sub3(p, a);
    float d1 = dot3(ab, ap), d2 = dot3(ac, ap);
    if (d1 <= 0.0f && d2 <= 0.0f) return dot3(ap, ap);
    V3 bp = sub3(p, b);
    float d3 = dot3(ab, bp), d4 = dot3(ac, bp);
    if (d3 >= 0.0f && d4 <= d3) return dot3(bp, bp);
    float vc = d1*d4 - d3*d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) return d2_seg(p, a, b);
    V3 cp = sub3(p, c);
    float d5 = dot3(ab, cp), d6 = dot3(ac, cp);
    if (d6 >= 0.0f && d5 <= d6) return dot3(cp, cp);
    float vb = d5*d2 - d1*d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) return d2_seg(p, a, c);
    float va = d3*d6 - d5*d4;
    if (va <= 0.0f && (d4-d3) >= 0.0f && (d5-d6) >= 0.0f) return d2_seg(p, b, c);
    V3 n = cross3(ab, ac);
    float nn = dot3(n, n);
    if (nn <= 0.0f) {              /* degenerate triangle: fall back to edges */
        float e0 = d2_seg(p, a, b), e1 = d2_seg(p, a, c), e2 = d2_seg(p, b, c);
        float m = e0 < e1 ? e0 : e1;
        return m < e2 ? m : e2;
    }
    float dist = dot3(n, ap);      /* inside face region: distance to plane */
    return dist*dist / nn;
}

/* characteristic vertex spacing from density; -1 when no measurable surface */
static float char_spacing(const Mesh* m, double* area_out) {
    double area = 0.0;
    for (int i = 0; i < m->nf; i++) {
        V3 a = m->v[m->f[i].a], b = m->v[m->f[i].b], c = m->v[m->f[i].c];
        if (!finite3(a) || !finite3(b) || !finite3(c)) continue;
        area += 0.5 * len3(cross3(sub3(b, a), sub3(c, a)));
    }
    if (area_out) *area_out = area;
    int nv = finite_verts(m);
    if (area <= 0.0 || nv == 0) return -1.0f;
    /* F ~= 2V (Euler), mean tri area A/2V, equilateral edge sqrt(4*area/sqrt(3)) */
    return (float)sqrt(2.0 * area / (1.7320508 * nv));
}

/* class 4: farthest unreferenced vertex, in units of s; 0 if none; -1 n/a */
static float iso_metric(const Mesh* m, float k, int* iso_count) {
    float s = char_spacing(m, NULL);
    *iso_count = 0;
    if (s <= 0.0f) return -1.0f;
    char ref[MAXV] = { 0 };
    for (int i = 0; i < m->nf; i++) {
        ref[m->f[i].a] = 1; ref[m->f[i].b] = 1; ref[m->f[i].c] = 1;
    }
    float maxd = 0.0f;
    for (int i = 0; i < m->nv; i++) {
        if (ref[i] || !finite3(m->v[i])) continue;
        float best = 1e30f;
        for (int j = 0; j < m->nf; j++) {
            V3 a = m->v[m->f[j].a], b = m->v[m->f[j].b], c = m->v[m->f[j].c];
            if (!finite3(a) || !finite3(b) || !finite3(c)) continue;
            float d2 = d2_tri(m->v[i], a, b, c);
            if (d2 < best) best = d2;
        }
        if (best >= 1e30f) continue;
        float ds = sqrtf(best) / s;
        if (ds > maxd) maxd = ds;
        if (ds > k) (*iso_count)++;
    }
    return maxd;
}

/* ---------------- the import-time screen ---------------- */

enum { SC_ZERO_VERTS = 1, SC_ZERO_VOLUME = 2, SC_NORMALS_ALIGNED = 4, SC_ISOLATED_VERTS = 8 };

#define SC_MIN_VERTS  3        /* fewer than 3 vertices cannot form even one face */
#define SC_BBOX_TAU   1e-4f    /* min extent / diagonal */
#define SC_NORMAL_TAU 0.95f    /* mean resultant length */
#define SC_ISO_K      0.5f     /* isolation distance = SC_ISO_K * s */

typedef struct { float bbox_tau, normal_tau, iso_k; } ScParams;

typedef struct {
    int flags;
    int finite;
    float bbox_ratio;   /* -1 n/a */
    float normal_r;     /* -1 n/a */
    float iso_ds;       /* -1 n/a */
    int iso_count;
} ScReport;

static void sc_run(const Mesh* m, const ScParams* p, ScReport* r) {
    memset(r, 0, sizeof(*r));
    r->bbox_ratio = -1.0f; r->normal_r = -1.0f; r->iso_ds = -1.0f;
    r->finite = finite_verts(m);
    if (r->finite < SC_MIN_VERTS) { r->flags |= SC_ZERO_VERTS; return; }
    r->bbox_ratio = bbox_ratio(m);
    if (r->bbox_ratio < p->bbox_tau) r->flags |= SC_ZERO_VOLUME;
    r->normal_r = normal_r(m);
    if (r->normal_r >= 0.0f && r->normal_r > p->normal_tau) r->flags |= SC_NORMALS_ALIGNED;
    if (m->nf == 0) { r->flags |= SC_ISOLATED_VERTS; return; }  /* no surface: every vertex isolated */
    r->iso_ds = iso_metric(m, p->iso_k, &r->iso_count);
    if (r->iso_count > 0) r->flags |= SC_ISOLATED_VERTS;
}

static const char* flags_str(int f) {
    static char b[80];
    b[0] = 0;
    if (!f) return "none";
    if (f & SC_ZERO_VERTS) strcat(b, "zero_verts ");
    if (f & SC_ZERO_VOLUME) strcat(b, "zero_volume ");
    if (f & SC_NORMALS_ALIGNED) strcat(b, "normals_aligned ");
    if (f & SC_ISOLATED_VERTS) strcat(b, "isolated_verts ");
    return b;
}

/* ---------------- test corpus: generators ---------------- */

static void compute_normals(Mesh* m) {
    for (int i = 0; i < m->nv; i++) m->n[i].x = m->n[i].y = m->n[i].z = 0.0f;
    for (int i = 0; i < m->nf; i++) {
        V3 a = m->v[m->f[i].a], b = m->v[m->f[i].b], c = m->v[m->f[i].c];
        if (!finite3(a) || !finite3(b) || !finite3(c)) continue;
        V3 n = cross3(sub3(b, a), sub3(c, a));   /* area-weighted */
        int idx[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
        for (int k = 0; k < 3; k++) {
            m->n[idx[k]].x += n.x; m->n[idx[k]].y += n.y; m->n[idx[k]].z += n.z;
        }
    }
    for (int i = 0; i < m->nv; i++) {
        float l = len3(m->n[i]);
        if (l > 1e-12f) { m->n[i].x /= l; m->n[i].y /= l; m->n[i].z /= l; }
    }
    m->has_normals = 1;
}

static void gen_box(Mesh* m, float sx, float sy, float sz) {
    static const float c[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    static const int t[12][3] = {
        {0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},
        {2,3,7},{2,7,6},{0,7,3},{0,4,7},{1,2,6},{1,6,5}
    };
    mesh_init(m);
    for (int i = 0; i < 8; i++) add_v(m, c[i][0]*sx, c[i][1]*sy, c[i][2]*sz);
    for (int i = 0; i < 12; i++) add_f(m, t[i][0], t[i][1], t[i][2]);
}

/* lat-long sphere (or cap when phimax < PI); poles duplicated per segment */
static void gen_uv_sphere(Mesh* m, float r, float phimax, int nr, int nseg) {
    mesh_init(m);
    for (int i = 0; i <= nr; i++) {
        float phi = phimax * i / nr;
        for (int j = 0; j < nseg; j++) {
            float th = 2.0f * PI * j / nseg;
            add_v(m, r*sinf(phi)*cosf(th), r*sinf(phi)*sinf(th), r*cosf(phi));
        }
    }
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nseg; j++) {
            int a = i*nseg + j, b = i*nseg + (j+1)%nseg;
            int c = a + nseg, d = b + nseg;
            add_f(m, a, c, b); add_f(m, b, c, d);
        }
}

static void gen_torus(Mesh* m, float R, float r, int nu, int nv) {
    mesh_init(m);
    for (int i = 0; i < nu; i++) {
        float u = 2.0f * PI * i / nu;
        for (int j = 0; j < nv; j++) {
            float w = 2.0f * PI * j / nv;
            add_v(m, (R + r*cosf(w))*cosf(u), (R + r*cosf(w))*sinf(u), r*sinf(w));
        }
    }
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nv; j++) {
            int a = i*nv + j, b = ((i+1)%nu)*nv + j;
            int c = i*nv + (j+1)%nv, d = ((i+1)%nu)*nv + (j+1)%nv;
            add_f(m, a, b, d); add_f(m, a, d, c);
        }
}

static unsigned rng_state = 12345u;
static float frand(void) {   /* deterministic [0,1) */
    rng_state = rng_state * 1664525u + 1013904223u;
    return (float)(rng_state >> 8) / 16777216.0f;
}

/* ---------------- test harness ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

typedef float (*Metric)(const Mesh*);
static float m_finite(const Mesh* m) { return (float)finite_verts(m); }
static float m_bbox(const Mesh* m) { return bbox_ratio(m); }
static float m_normr(const Mesh* m) { return normal_r(m); }
static float m_iso(const Mesh* m) { int c; return iso_metric(m, SC_ISO_K, &c); }

/* print FP/MISS vs threshold; the zero-error plateau IS the tradeoff picture */
static void sweep(const char* tag, Metric fn, const float* taus, int nt, int less,
                  float chosen, const Mesh* L, int nl, const Mesh* B, int nb) {
    printf("== sweep %s (flag if metric %s tau) ==\n", tag, less ? "<" : ">");
    int lo = -1, hi = -1, chosen_ok = 0;
    for (int t = 0; t < nt; t++) {
        float tau = taus[t];
        int fp = 0, miss = 0;
        for (int i = 0; i < nl; i++) {
            float v = fn(&L[i]);
            if (v >= 0.0f && (less ? v < tau : v > tau)) fp++;
        }
        for (int i = 0; i < nb; i++) {
            float v = fn(&B[i]);
            if (v >= 0.0f && !(less ? v < tau : v > tau)) miss++;
        }
        int ok = (fp == 0 && miss == 0);
        if (ok) { if (lo < 0) lo = t; hi = t; if (tau == chosen) chosen_ok = 1; }
        printf("  tau=%-10.4g FP=%d/%d MISS=%d/%d%s\n",
               tau, fp, nl, miss, nb, tau == chosen ? "  <-- chosen" : "");
    }
    if (lo >= 0) printf("  zero-error plateau: [%g, %g]\n", taus[lo], taus[hi]);
    char nm[96];
    snprintf(nm, sizeof nm, "sweep_%s_chosen_in_plateau", tag);
    check(nm, chosen_ok);
}

int main(void) {
    const ScParams P = { SC_BBOX_TAU, SC_NORMAL_TAU, SC_ISO_K };

    /* ---- legitimate batch: must produce zero false positives ---- */
    Mesh L[6];
    const char* ln[6];

    gen_box(&L[0], 1, 1, 1); compute_normals(&L[0]);
    ln[0] = "cube";

    gen_uv_sphere(&L[1], 1, PI, 16, 16); compute_normals(&L[1]);
    ln[1] = "sphere";

    gen_box(&L[2], 1, 1, 3e-3f); compute_normals(&L[2]);   /* sheet metal ~2e-3 ratio */
    ln[2] = "plate_3e-3";

    gen_uv_sphere(&L[3], 2, PI/3, 12, 16); compute_normals(&L[3]);  /* 60 deg cap: R~0.8 */
    ln[3] = "cap_60deg";

    gen_uv_sphere(&L[4], 1, PI, 16, 16);
    {   /* leftover vertices from editing: on-surface (d~sagitta) and 0.08 s off */
        float s = char_spacing(&L[4], NULL);
        add_v(&L[4], sinf(1.0f)*cosf(2.0f), sinf(1.0f)*sinf(2.0f), cosf(1.0f));
        add_v(&L[4], sinf(2.0f)*cosf(4.0f), sinf(2.0f)*sinf(4.0f), cosf(2.0f));
        add_v(&L[4], sinf(2.5f)*cosf(1.0f), sinf(2.5f)*sinf(1.0f), cosf(2.5f));
        float d = 1.0f + 0.08f * s;
        add_v(&L[4], d*sinf(0.7f), 0, d*cosf(0.7f));
        add_v(&L[4], 0, d*sinf(2.8f), d*cosf(2.8f));
    }
    compute_normals(&L[4]);
    ln[4] = "sphere_stray";

    gen_torus(&L[5], 1, 0.4f, 16, 8); compute_normals(&L[5]);
    ln[5] = "torus";

    printf("-- legit batch (expect flags=none) --\n");
    for (int i = 0; i < 6; i++) {
        ScReport r;
        sc_run(&L[i], &P, &r);
        printf("%-14s finite=%3d bbox=%.3e R=%.3f iso_ds=%.2f flags=%s\n",
               ln[i], r.finite, r.bbox_ratio, r.normal_r, r.iso_ds, flags_str(r.flags));
        char nm[96];
        snprintf(nm, sizeof nm, "legit_%s_clean", ln[i]);
        check(nm, r.flags == 0);
    }

    /* ---- injected batch: two manifestations per class ---- */
    Mesh B[8];
    const char* bn[8];
    int bexp[8];

    mesh_init(&B[0]);                                    /* 1a: empty buffer */
    bn[0] = "empty_buffer";    bexp[0] = SC_ZERO_VERTS;

    gen_box(&B[1], 1, 1, 1);                             /* 1b: truncated/garbled stream: */
    for (int i = 2; i < B[1].nv; i++)                    /*     2 finite verts, rest NaN */
        B[1].v[i].x = B[1].v[i].y = B[1].v[i].z = NAN;
    bn[1] = "truncated_nan";   bexp[1] = SC_ZERO_VERTS;

    gen_box(&B[2], 0, 0, 0);                             /* 2a: collapsed to a point */
    bn[2] = "point_collapse";  bexp[2] = SC_ZERO_VOLUME;

    gen_box(&B[3], 1, 1, 2e-5f);                         /* 2b: flattened to ~1e-5 (near-zero) */
    bn[3] = "near_flat";       bexp[3] = SC_ZERO_VOLUME;

    gen_box(&B[4], 1, 1, 1);                             /* 3a: normal buffer filled with */
    B[4].has_normals = 1;                                /*     one constant */
    for (int i = 0; i < B[4].nv; i++) B[4].n[i] = (V3){ 0, 0, 1 };
    bn[4] = "normals_const";   bexp[4] = SC_NORMALS_ALIGNED;

    gen_box(&B[5], 1, 1, 1);                             /* 3b: constant + 1% jitter */
    B[5].has_normals = 1;
    for (int i = 0; i < B[5].nv; i++) {
        V3 nn = { 0.01f*(2*frand()-1), 0.01f*(2*frand()-1), 1.0f };
        float l = len3(nn);
        B[5].n[i] = (V3){ nn.x/l, nn.y/l, nn.z/l };
    }
    bn[5] = "normals_jitter";  bexp[5] = SC_NORMALS_ALIGNED;

    gen_box(&B[6], 1, 1, 1);                             /* 4a: garbage vertices far away */
    for (int j = 0; j < 5; j++)
        add_v(&B[6], 10.0f + j, 10.0f - 0.3f*j, 10.0f + 0.2f*j);
    compute_normals(&B[6]);
    bn[6] = "far_garbage";     bexp[6] = SC_ISOLATED_VERTS;

    gen_box(&B[7], 1, 1, 1);                             /* 4b: strays one lattice step off (3 s) */
    {
        float s = char_spacing(&B[7], NULL);
        add_v(&B[7], 0.5f, 0.5f, 1.0f + 3.0f*s);
        add_v(&B[7], 0.5f, 0.5f, -3.0f*s);
        add_v(&B[7], 0.5f, 1.0f + 3.0f*s, 0.5f);
        add_v(&B[7], -3.0f*s, 0.5f, 0.5f);
    }
    compute_normals(&B[7]);
    bn[7] = "near_stray";      bexp[7] = SC_ISOLATED_VERTS;

    printf("-- injected batch (expect its class bit, and only it) --\n");
    for (int i = 0; i < 8; i++) {
        ScReport r;
        sc_run(&B[i], &P, &r);
        printf("%-14s finite=%3d bbox=%.3e R=%.3f iso_ds=%.2f flags=%s\n",
               bn[i], r.finite, r.bbox_ratio, r.normal_r, r.iso_ds, flags_str(r.flags));
        char nm[96];
        snprintf(nm, sizeof nm, "bad_%s_caught", bn[i]);
        check(nm, (r.flags & bexp[i]) == bexp[i]);
        snprintf(nm, sizeof nm, "bad_%s_exact", bn[i]);
        check(nm, r.flags == bexp[i]);
    }

    /* ---- threshold sweeps: the FP/FN tradeoff for each knob ---- */
    static const float k_tau[]    = { 1, 2, 3, 4, 8, 16 };
    static const float bbox_tau[] = { 1e-7f, 1e-6f, 1e-5f, 1e-4f, 1e-3f, 1e-2f, 1e-1f };
    static const float r_tau[]    = { 0.5f, 0.6f, 0.7f, 0.75f, 0.85f, 0.9f, 0.95f, 0.99f, 0.999f, 0.99999f };
    static const float iso_k[]    = { 0.02f, 0.05f, 0.1f, 0.2f, 0.5f, 1.0f, 2.0f, 4.0f };

    sweep("zero_verts", m_finite, k_tau, 6, 1, SC_MIN_VERTS, L, 6, &B[0], 2);
    sweep("zero_volume", m_bbox, bbox_tau, 7, 1, SC_BBOX_TAU, L, 6, &B[2], 2);
    sweep("normals", m_normr, r_tau, 10, 0, SC_NORMAL_TAU, L, 6, &B[4], 2);
    sweep("isolated", m_iso, iso_k, 8, 0, SC_ISO_K, L, 6, &B[6], 2);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
