
/* lod: skeleton-locked LOD chains.
 *
 * Hard rule: every vertex of level k must sit on an original vertex or edge
 * of level k-1.  Placement is classified geometrically as vertex/edge/face.
 * Second criterion: sampled symmetric Hausdorff ("max surface deviation")
 * between adjacent levels must stay under the model target eps.
 * A deviation-only decimator (vertices free to drift, UV left on the stale
 * lattice parameter) is run at the same sizes as a diagnostic baseline.
 *
 * `make test` builds and runs the self-checking report; `./lod --csv` emits
 * the 18 per-level metric rows.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TWOPI 6.28318530717958647692

#define FINE_N 65
#define MAXN   65
#define MAXV   (FINE_N * FINE_N)
#define MAXF   ((FINE_N - 1) * (FINE_N - 1) * 2)

#define TEX_SIZE 256.0
#define LOC_TOL  1e-6      /* world units; locked refs match to ~1e-15 */
#define DEV_K    64        /* DEV_K^2 domain samples per side */
#define UV_K     96

typedef struct { double x, y, z; } V3;
typedef struct { double u, v; } UV2;
typedef struct { int a, b, c; } Tri;

typedef struct {
    int nx, ny;
    double X[MAXN], Y[MAXN];
    V3  p[MAXV];
    UV2 uv[MAXV];
    Tri f[MAXF];
    int nv, nf;
} Mesh;

enum { LOC_VERTEX = 0, LOC_EDGE = 1, LOC_FACE = 2 };

static int failures = 0;
static void check(const char* name, int ok) {
    printf("check %s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---------------- models: 3 surfaces, 3 deviation targets ---------------- */

static const char* MODEL_NAME[3] = { "M1_dome", "M2_wave", "M3_ridge" };
static const double EPS_TARGET[3] = { 0.0050, 0.0200, 0.0500 };
static const double UV_AMP[3]     = { 0.030, 0.050, 0.075 };

/* one density ladder per target tier: tight -> dense, loose -> coarse */
static const int LADDER[3][3][2] = {
    { {43,33}, {25,19}, {15,11} },
    { {35,27}, {21,16}, {12, 9} },
    { {25,19}, {15,11}, { 9, 7} }
};

static double surface(int m, double x, double y) {
    double xx = 2.0 * x - 1.0, yy = 2.0 * y - 1.0;
    if (m == 0) return 0.18 * (1.0 - xx * xx - yy * yy) + 0.05;
    if (m == 1) return 0.10 * sin(TWOPI * x) * cos(TWOPI * y);
    return 0.07 * sin(4.0 * PI * x + 0.6 * sin(TWOPI * y)) * cos(TWOPI * y);
}

static void uvmap(int m, double x, double y, UV2* o) {
    double a = UV_AMP[m];
    o->u = x + a * sin(TWOPI * y);
    o->v = y + a * sin(TWOPI * x + 0.7);
}

/* ---------------- structured heightfield grid ---------------- */

static int grid_tri_count(int nx, int ny) { return 2 * (nx - 1) * (ny - 1); }

static void build_fine(Mesh* M, int model) {
    memset(M, 0, sizeof(*M));
    M->nx = M->ny = FINE_N;
    M->nv = FINE_N * FINE_N;
    M->nf = grid_tri_count(FINE_N, FINE_N);
    for (int i = 0; i < FINE_N; i++) M->X[i] = (double)i / (FINE_N - 1);
    for (int j = 0; j < FINE_N; j++) M->Y[j] = (double)j / (FINE_N - 1);
    for (int j = 0; j < FINE_N; j++) {
        for (int i = 0; i < FINE_N; i++) {
            int k = j * FINE_N + i;
            double x = M->X[i], y = M->Y[j];
            M->p[k].x = x; M->p[k].y = y; M->p[k].z = surface(model, x, y);
            uvmap(model, x, y, &M->uv[k]);
        }
    }
    int n = 0;
    for (int j = 0; j < FINE_N - 1; j++) {
        for (int i = 0; i < FINE_N - 1; i++) {
            int a = j * FINE_N + i, b = a + 1, c = a + FINE_N, d = c + 1;
            M->f[n++] = (Tri){ a, c, b };   /* lower-left triangle  */
            M->f[n++] = (Tri){ b, c, d };   /* upper-right triangle */
        }
    }
}

/* locate (x,y) in the grid; barycentric weights match the diagonal above */
static void grid_weights(const Mesh* M, double x, double y,
                         int id[3], double w[3]) {
    int i = 0, j = 0;
    while (i < M->nx - 2 && x > M->X[i + 1]) i++;
    while (j < M->ny - 2 && y > M->Y[j + 1]) j++;
    if (i > M->nx - 2) i = M->nx - 2;
    if (j > M->ny - 2) j = M->ny - 2;
    double dx = M->X[i + 1] - M->X[i], dy = M->Y[j + 1] - M->Y[j];
    double fx = (x - M->X[i]) / dx, fy = (y - M->Y[j]) / dy;
    if (fx < 0) fx = 0;
    if (fx > 1) fx = 1;
    if (fy < 0) fy = 0;
    if (fy > 1) fy = 1;
    int i00 = j * M->nx + i, i10 = i00 + 1;
    int i01 = i00 + M->nx,     i11 = i01 + 1;
    if (fx + fy <= 1.0) {
        id[0] = i00; id[1] = i01; id[2] = i10;
        w[0] = 1.0 - fx - fy; w[1] = fy; w[2] = fx;
    } else {
        id[0] = i10; id[1] = i01; id[2] = i11;
        w[0] = 1.0 - fy; w[1] = 1.0 - fx; w[2] = fx + fy - 1.0;
    }
}

static V3 grid_pos(const Mesh* M, double x, double y) {
    int id[3]; double w[3];
    grid_weights(M, x, y, id, w);
    V3 r = { 0,0,0 };
    for (int k = 0; k < 3; k++) {
        r.x += w[k] * M->p[id[k]].x;
        r.y += w[k] * M->p[id[k]].y;
        r.z += w[k] * M->p[id[k]].z;
    }
    return r;
}

static UV2 grid_uv(const Mesh* M, double x, double y) {
    int id[3]; double w[3];
    grid_weights(M, x, y, id, w);
    UV2 r = { 0,0 };
    for (int k = 0; k < 3; k++) {
        r.u += w[k] * M->uv[id[k]].u;
        r.v += w[k] * M->uv[id[k]].v;
    }
    return r;
}

/* ---------------- locked decimator: refs are parent vertices/edges ---------------- */

typedef struct { int j; int edge; double t; } XRef;

static V3 lerp3(const V3* a, const V3* b, double t) {
    return (V3){ a->x + (b->x - a->x) * t,
                 a->y + (b->y - a->y) * t,
                 a->z + (b->z - a->z) * t };
}
static UV2 lerp2(const UV2* a, const UV2* b, double t) {
    return (UV2){ a->u + (b->u - a->u) * t,
                  a->v + (b->v - a->v) * t };
}

static void coarsen_locked(const Mesh* P, Mesh* C, int cx, int cy) {
    memset(C, 0, sizeof(*C));
    C->nx = cx; C->ny = cy;
    XRef xr[MAXN];
    int yidx[MAXN];

    for (int i = 0; i < cx; i++) {
        double q = (double)i * (P->nx - 1) / (cx - 1);
        int j = (int)floor(q);
        double f = q - j;
        if (j >= P->nx - 1) { j = P->nx - 2; f = 1.0; }
        if (f < 0.18) {
            xr[i] = (XRef){ j, 0, 0.0 };            /* snap to vertex j   */
        } else if (f > 0.82) {
            int jj = j + 1;
            if (jj > P->nx - 1) jj = P->nx - 1;
            xr[i] = (XRef){ jj, 0, 0.0 };            /* snap to vertex j+1 */
        } else {
            xr[i] = (XRef){ j, 1, f };               /* interior of edge j */
        }
    }
    for (int r = 0; r < cy; r++) {
        double q = (double)r * (P->ny - 1) / (cy - 1);
        int k = (int)lround(q);
        if (k < 0) k = 0;
        if (k > P->ny - 1) k = P->ny - 1;
        yidx[r] = k;                                 /* y always a vertex row */
    }

    C->nv = cx * cy;
    C->nf = grid_tri_count(cx, cy);
    for (int i = 0; i < cx; i++) {
        int j = xr[i].j;
        double t = xr[i].edge ? xr[i].t : 0.0;
        int src = yidx[0] * P->nx + j;
        V3 xp = xr[i].edge
            ? lerp3(&P->p[src], &P->p[src + 1], t)
            : P->p[src];
        C->X[i] = xp.x;
    }
    for (int r = 0; r < cy; r++) C->Y[r] = P->Y[yidx[r]];

    for (int r = 0; r < cy; r++) {
        int k = yidx[r];
        for (int i = 0; i < cx; i++) {
            int o = r * cx + i;
            int j = xr[i].j;
            if (!xr[i].edge) {
                int s = k * P->nx + j;
                C->p[o] = P->p[s];
                C->uv[o] = P->uv[s];
            } else {
                int a = k * P->nx + j, b = a + 1;
                C->p[o]  = lerp3(&P->p[a], &P->p[b], xr[i].t);
                C->uv[o] = lerp2(&P->uv[a], &P->uv[b], xr[i].t);
            }
        }
    }
    int n = 0;
    for (int r = 0; r < cy - 1; r++) {
        for (int i = 0; i < cx - 1; i++) {
            int a = r * cx + i, b = a + 1, cc = a + cx, d = cc + 1;
            C->f[n++] = (Tri){ a, cc, b };
            C->f[n++] = (Tri){ b, cc, d };
        }
    }
}

/* ---------------- deviation-only decimator: free tangential drift, stale UV ---------------- */

static void coarsen_free(const Mesh* P, Mesh* C, int cx, int cy,
                         int model, double drift) {
    memset(C, 0, sizeof(*C));
    C->nx = cx; C->ny = cy;
    double xs[MAXN], x0s[MAXN], ys[MAXN], y0s[MAXN];

    for (int i = 0; i < cx; i++) {
        double q = (double)i * (P->nx - 1) / (cx - 1);
        int j = (int)floor(q);
        double f = q - j;
        if (j >= P->nx - 1) { j = P->nx - 2; f = 1.0; }
        double x0 = P->X[j] + f * (P->X[j + 1] - P->X[j]);
        double hx = P->X[j + 1] - P->X[j];
        int boundary = (i == 0 || i == cx - 1);
        double dx = boundary ? 0.0
                             : drift * hx * sin(1.7 * i + 0.9 * cy);
        x0s[i] = x0; xs[i] = x0 + dx;
    }
    for (int r = 0; r < cy; r++) {
        double q = (double)r * (P->ny - 1) / (cy - 1);
        int j = (int)floor(q);
        double f = q - j;
        if (j >= P->ny - 1) { j = P->ny - 2; f = 1.0; }
        double y0 = P->Y[j] + f * (P->Y[j + 1] - P->Y[j]);
        double hy = P->Y[j + 1] - P->Y[j];
        int boundary = (r == 0 || r == cy - 1);
        double dy = boundary ? 0.0
                             : drift * hy * cos(1.3 * r + 1.1 * cx);
        y0s[r] = y0; ys[r] = y0 + dy;
    }
    for (int i = 0; i < cx; i++) C->X[i] = xs[i];
    for (int r = 0; r < cy; r++) C->Y[r] = ys[r];

    C->nv = cx * cy;
    C->nf = grid_tri_count(cx, cy);
    for (int r = 0; r < cy; r++) {
        for (int i = 0; i < cx; i++) {
            int o = r * cx + i;
            C->p[o] = (V3){ xs[i], ys[r], surface(model, xs[i], ys[r]) };
            /* UV channel stays bound to the pre-drift lattice parameter:
             * this is the correspondence bug a deviation-only pass leaves. */
            uvmap(model, x0s[i], y0s[r], &C->uv[o]);
        }
    }
    int n = 0;
    for (int r = 0; r < cy - 1; r++) {
        for (int i = 0; i < cx - 1; i++) {
            int a = r * cx + i, b = a + 1, cc = a + cx, d = cc + 1;
            C->f[n++] = (Tri){ a, cc, b };
            C->f[n++] = (Tri){ b, cc, d };
        }
    }
}

/* ---------------- placement classification vs parent skeleton ---------------- */

static double d3sq(const V3* a, const V3* b) {
    double x = a->x - b->x, y = a->y - b->y, z = a->z - b->z;
    return x * x + y * y + z * z;
}

static double seg_dsq(const V3* p, const V3* a, const V3* b) {
    V3 ab = { b->x - a->x, b->y - a->y, b->z - a->z };
    double l2 = ab.x * ab.x + ab.y * ab.y + ab.z * ab.z;
    double t = 0.0;
    if (l2 > 0.0) {
        t = ((p->x - a->x) * ab.x + (p->y - a->y) * ab.y
             + (p->z - a->z) * ab.z) / l2;
        if (t < 0) t = 0; else if (t > 1) t = 1;
    }
    V3 q = { a->x + t * ab.x, a->y + t * ab.y, a->z + t * ab.z };
    return d3sq(p, &q);
}

static void classify_placement(const Mesh* C, const Mesh* P, int cnt[3]) {
    cnt[0] = cnt[1] = cnt[2] = 0;
    double tol2 = LOC_TOL * LOC_TOL;
    for (int k = 0; k < C->nv; k++) {
        const V3* p = &C->p[k];
        double dv2 = 1e300;
        for (int j = 0; j < P->nv; j++) {
            double q = d3sq(p, &P->p[j]);
            if (q < dv2) dv2 = q;
        }
        if (dv2 <= tol2) { cnt[LOC_VERTEX]++; continue; }

        double de2 = 1e300;
        for (int ti = 0; ti < P->nf; ti++) {
            int ix[3] = { P->f[ti].a, P->f[ti].b, P->f[ti].c };
            for (int e = 0; e < 3; e++) {
                double q = seg_dsq(p, &P->p[ix[e]], &P->p[ix[(e + 1) % 3]]);
                if (q < de2) de2 = q;
            }
        }
        if (de2 <= tol2) cnt[LOC_EDGE]++;
        else             cnt[LOC_FACE]++;
    }
}

/* ---------------- exact point/triangle distance (Ericson) ---------------- */

static V3 closest_pt_tri(const V3* p, const V3* a, const V3* b, const V3* c) {
    V3 ab = { b->x-a->x, b->y-a->y, b->z-a->z };
    V3 ac = { c->x-a->x, c->y-a->y, c->z-a->z };
    V3 ap = { p->x-a->x, p->y-a->y, p->z-a->z };
    double d1 = ab.x*ap.x + ab.y*ap.y + ab.z*ap.z;
    double d2 = ac.x*ap.x + ac.y*ap.y + ac.z*ap.z;
    if (d1 <= 0 && d2 <= 0) return *a;

    V3 bp = { p->x-b->x, p->y-b->y, p->z-b->z };
    double d3 = ab.x*bp.x + ab.y*bp.y + ab.z*bp.z;
    double d4 = ac.x*bp.x + ac.y*bp.y + ac.z*bp.z;
    if (d3 >= 0 && d4 <= d3) return *b;

    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        return (V3){ a->x + v*ab.x, a->y + v*ab.y, a->z + v*ab.z };
    }
    V3 cp = { p->x-c->x, p->y-c->y, p->z-c->z };
    double d5 = ab.x*cp.x + ab.y*cp.y + ab.z*cp.z;
    double d6 = ac.x*cp.x + ac.y*cp.y + ac.z*cp.z;
    if (d6 >= 0 && d5 <= d6) return *c;

    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        return (V3){ a->x + w*ac.x, a->y + w*ac.y, a->z + w*ac.z };
    }
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        V3 bc = { c->x-b->x, c->y-b->y, c->z-b->z };
        return (V3){ b->x + w*bc.x, b->y + w*bc.y, b->z + w*bc.z };
    }
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom, w = vc * denom;
    return (V3){ a->x + ab.x*v + ac.x*w,
                 a->y + ab.y*v + ac.y*w,
                 a->z + ab.z*v + ac.z*w };
}

static double pt_mesh_dsq(const V3* p, const Mesh* T) {
    double best = 1e300;
    for (int ti = 0; ti < T->nf; ti++) {
        const Tri* tr = &T->f[ti];
        V3 q = closest_pt_tri(p, &T->p[tr->a], &T->p[tr->b], &T->p[tr->c]);
        double dd = d3sq(p, &q);
        if (dd < best) best = dd;
    }
    return best;
}

/* sampled symmetric Hausdorff: domain grid samples + every vertex */
static double dev_one_side(const Mesh* S, const Mesh* T, int K) {
    double worst = 0.0;
    for (int j = 0; j <= K; j++) {
        double y = (double)j / K;
        for (int i = 0; i <= K; i++) {
            double x = (double)i / K;
            V3 p = grid_pos(S, x, y);
            double dd = pt_mesh_dsq(&p, T);
            if (dd > worst) worst = dd;
        }
    }
    for (int k = 0; k < S->nv; k++) {
        double dd = pt_mesh_dsq(&S->p[k], T);
        if (dd > worst) worst = dd;
    }
    return worst;
}

static double surface_dev(const Mesh* A, const Mesh* B) {
    double d1 = dev_one_side(A, B, DEV_K);
    double d2 = dev_one_side(B, A, DEV_K);
    return sqrt(d1 > d2 ? d1 : d2);
}

/* ---------------- UV / texture-shift metrics ---------------- */

static double uv_field_shift(const Mesh* A, const Mesh* B, int K,
                             double* rms_px) {
    double worst = 0.0, sum = 0.0;
    long n = 0;
    for (int j = 0; j <= K; j++) {
        double y = (double)j / K;
        for (int i = 0; i <= K; i++) {
            double x = (double)i / K;
            UV2 a = grid_uv(A, x, y), b = grid_uv(B, x, y);
            double du = (a.u - b.u) * TEX_SIZE;
            double dv = (a.v - b.v) * TEX_SIZE;
            double dd = sqrt(du * du + dv * dv);
            if (dd > worst) worst = dd;
            sum += dd * dd; n++;
        }
    }
    *rms_px = sqrt(sum / n);
    return worst;
}

/* mismatch at the new level's own vertices vs the parent UV field */
static double uv_anchor_shift(const Mesh* Child, const Mesh* Parent) {
    double worst = 0.0;
    for (int k = 0; k < Child->nv; k++) {
        UV2 pu = grid_uv(Parent, Child->p[k].x, Child->p[k].y);
        double du = (Child->uv[k].u - pu.u) * TEX_SIZE;
        double dv = (Child->uv[k].v - pu.v) * TEX_SIZE;
        double dd = sqrt(du * du + dv * dv);
        if (dd > worst) worst = dd;
    }
    return worst;
}

/* ---------------- experiment ---------------- */

typedef struct {
    int nv;
    int cnt[3];
    double pct[3];
    double dmax, uvmax, uvrms, anchor;
    int dev_pass;
} LevelMet;

typedef struct {
    double eps;
    LevelMet lk[3], fr[3];
    double dmax_lk, uvmax_lk, anchor_lk;
    double dmax_fr, uvmax_fr, anchor_fr;
    int pool_lk[3], pool_fr[3], pool_nv_lk, pool_nv_fr;
    int hard_pass, qualified, free_dev_pass;
} ModelMet;

static Mesh g_fine, g_lk[3], g_fr[3];
static ModelMet MM[3];

static LevelMet measure_level(const Mesh* Child, const Mesh* Parent,
                              double eps) {
    LevelMet L;
    memset(&L, 0, sizeof(L));
    L.nv = Child->nv;
    classify_placement(Child, Parent, L.cnt);
    for (int k = 0; k < 3; k++) L.pct[k] = 100.0 * L.cnt[k] / Child->nv;
    L.dmax = surface_dev(Parent, Child);
    L.uvmax = uv_field_shift(Parent, Child, UV_K, &L.uvrms);
    L.anchor = uv_anchor_shift(Child, Parent);
    L.dev_pass = (L.dmax <= eps);
    return L;
}

static void run_model(int m) {
    ModelMet* R = &MM[m];
    memset(R, 0, sizeof(*R));
    R->eps = EPS_TARGET[m];
    build_fine(&g_fine, m);

    const Mesh* pl = &g_fine;
    const Mesh* pf = &g_fine;
    for (int l = 0; l < 3; l++) {
        int cx = LADDER[m][l][0], cy = LADDER[m][l][1];
        coarsen_locked(pl, &g_lk[l], cx, cy);
        R->lk[l] = measure_level(&g_lk[l], pl, R->eps);

        coarsen_free(pf, &g_fr[l], cx, cy, m, 0.42);
        R->fr[l] = measure_level(&g_fr[l], pf, R->eps);

        R->dmax_lk = fmax(R->dmax_lk, R->lk[l].dmax);
        R->uvmax_lk = fmax(R->uvmax_lk, R->lk[l].uvmax);
        R->anchor_lk = fmax(R->anchor_lk, R->lk[l].anchor);
        R->dmax_fr = fmax(R->dmax_fr, R->fr[l].dmax);
        R->uvmax_fr = fmax(R->uvmax_fr, R->fr[l].uvmax);
        R->anchor_fr = fmax(R->anchor_fr, R->fr[l].anchor);
        for (int k = 0; k < 3; k++) {
            R->pool_lk[k] += R->lk[l].cnt[k];
            R->pool_fr[k] += R->fr[l].cnt[k];
        }
        R->pool_nv_lk += R->lk[l].nv;
        R->pool_nv_fr += R->fr[l].nv;

        pl = &g_lk[l];
        pf = &g_fr[l];
    }
    R->hard_pass = 1;
    for (int l = 0; l < 3; l++)
        if (R->lk[l].cnt[LOC_FACE] != 0) R->hard_pass = 0;
    R->qualified = R->hard_pass && (R->dmax_lk <= R->eps);
    R->free_dev_pass = (R->dmax_fr <= R->eps);
}

/* ---------------- output ---------------- */

static void print_row(const char* tag, int lvl, const LevelMet* L) {
    printf("  L%d %-6s n=%4d  V=%5.1f E=%5.1f F=%5.1f  dmax=%.5f %s"
           "  uvmax=%6.2f rms=%5.2f anchor=%7.3f\n",
           lvl + 1, tag, L->nv, L->pct[0], L->pct[1], L->pct[2],
           L->dmax, L->dev_pass ? "PASS" : "over",
           L->uvmax, L->uvrms, L->anchor);
}

static void print_report(void) {
    printf("lod: skeleton-locked LOD chains; fine grid %dx%d (%d v, %d tri),"
           " tex=%dpx, dev samples=%d^2+verts\n",
           FINE_N, FINE_N, MAXV, MAXF, (int)TEX_SIZE, DEV_K);
    printf("targets eps = %.4f / %.4f / %.4f\n\n",
           EPS_TARGET[0], EPS_TARGET[1], EPS_TARGET[2]);

    for (int m = 0; m < 3; m++) {
        const ModelMet* R = &MM[m];
        printf("== %s  eps=%.4f  ladder", MODEL_NAME[m], R->eps);
        for (int l = 0; l < 3; l++)
            printf(" %dx%d", LADDER[m][l][0], LADDER[m][l][1]);
        printf(" ==\n");
        for (int l = 0; l < 3; l++) {
            print_row("lock", l, &R->lk[l]);
            print_row("free", l, &R->fr[l]);
        }
        double pv = 100.0 * R->pool_lk[0] / R->pool_nv_lk;
        double pe = 100.0 * R->pool_lk[1] / R->pool_nv_lk;
        double pf = 100.0 * R->pool_lk[2] / R->pool_nv_lk;
        printf("  pooled lock V/E/F = %.1f/%.1f/%.1f  chain dmax=%.5f"
               " uvmax=%.2f anchor=%.4f  qualified=%d\n\n",
               pv, pe, pf, R->dmax_lk, R->uvmax_lk, R->anchor_lk,
               R->qualified);
    }

    printf("TWELVE_ORDER: each model group = "
           "[vertex_pct_pooled, edge_pct_pooled, dmax_chain, uvmax_px]\n");
    printf("TWELVE: ");
    for (int m = 0; m < 3; m++) {
        const ModelMet* R = &MM[m];
        double pv = 100.0 * R->pool_lk[0] / R->pool_nv_lk;
        double pe = 100.0 * R->pool_lk[1] / R->pool_nv_lk;
        printf("%s%.1f,%.1f,%.5f,%.2f", m ? "; " : "",
               pv, pe, R->dmax_lk, R->uvmax_lk);
    }
    printf("\n");

    printf("TWELVE_SCORECARD_ORDER: [face_pct_peak, dmax_chain, uvmax_px,"
           " qualified_both]\n");
    printf("TWELVE_SCORECARD: ");
    for (int m = 0; m < 3; m++) {
        const ModelMet* R = &MM[m];
        double fpeak = 0;
        for (int l = 0; l < 3; l++)
            fpeak = fmax(fpeak, R->lk[l].pct[2]);
        printf("%s%.1f,%.5f,%.2f,%d", m ? "; " : "",
               fpeak, R->dmax_lk, R->uvmax_lk, R->qualified);
    }
    printf("\n");

    /* conditional diagnostic: M3 deviation passes but texture jumps */
    const ModelMet* D = &MM[2];
    printf("TWELVE_DIAG_ORDER: M3 three levels x [V%%,E%%,F%%,uvmax_px]\n");
    printf("TWELVE_DIAG: ");
    for (int l = 0; l < 3; l++) {
        const LevelMet* L = &D->fr[l];
        printf("%s%.1f,%.1f,%.1f,%.2f", l ? "; " : "",
               L->pct[0], L->pct[1], L->pct[2], L->uvmax);
    }
    printf("\n\n");

    /* ---- self checks ---- */
    check("fine_vertex_count", MAXV == 4225);
    check("fine_tri_count", g_fine.nf == 8192);

    /* classifier unit geometry: corner vertex, edge midpoint, face centroid */
    Mesh U; memset(&U, 0, sizeof(U));
    U.nx = U.ny = 2; U.nv = 4;
    U.p[0] = (V3){0,0,0}; U.p[1] = (V3){1,0,0};
    U.p[2] = (V3){0,1,0}; U.p[3] = (V3){1,1,0};
    U.f[0] = (Tri){0, 2, 1}; U.f[1] = (Tri){1, 2, 3}; U.nf = 2;
    Mesh Q; memset(&Q, 0, sizeof(Q)); Q.nv = 3;
    Q.p[0] = U.p[0];
    Q.p[1] = (V3){0.5, 0.0, 0.0};
    Q.p[2] = (V3){1.0 / 3.0, 1.0 / 3.0, 0.0};
    int uc[3]; classify_placement(&Q, &U, uc);
    check("classify_vertex", uc[LOC_VERTEX] >= 1);
    check("classify_edge", uc[LOC_EDGE] >= 1);
    check("classify_face", uc[LOC_FACE] >= 1);

    for (int m = 0; m < 3; m++) {
        char nm[64];
        snprintf(nm, sizeof(nm), "m%d_locked_zero_face", m + 1);
        check(nm, MM[m].hard_pass);
        snprintf(nm, sizeof(nm), "m%d_locked_dev_under_target", m + 1);
        check(nm, MM[m].dmax_lk <= MM[m].eps);
        snprintf(nm, sizeof(nm), "m%d_locked_anchor_exact", m + 1);
        check(nm, MM[m].anchor_lk < 1e-3);
        snprintf(nm, sizeof(nm), "m%d_chain_qualified", m + 1);
        check(nm, MM[m].qualified == 1);
    }

    check("all_locked_qualified",
          MM[0].qualified && MM[1].qualified && MM[2].qualified);

    /* deviation-only arm: geometric budget still met ... */
    check("diag_m3_free_dev_pass", MM[2].free_dev_pass);
    /* ... but vertices sit in faces and the texture jumps */
    double m3_fpeak = 0, m3_uv = 0;
    for (int l = 0; l < 3; l++) {
        m3_fpeak = fmax(m3_fpeak, MM[2].fr[l].pct[2]);
        m3_uv = fmax(m3_uv, MM[2].fr[l].uvmax);
    }
    check("diag_m3_free_face_heavy", m3_fpeak >= 50.0);
    check("diag_m3_free_uv_jump", m3_uv >= 8.0);
    check("diag_face_explains_jump",
          MM[2].free_dev_pass && m3_fpeak >= 50.0 && m3_uv >= 8.0);
    check("diag_free_worse_than_locked", MM[2].uvmax_fr > MM[2].uvmax_lk);

    printf("\nfailures=%d\n", failures);
    printf("all_ok=%d\n", failures == 0);
}

static void print_csv(void) {
    printf("model,eps,arm,level,nv,pct_v,pct_e,pct_f,dmax,dev_pass,"
           "uvmax_px,uvrms_px,anchor_px\n");
    for (int m = 0; m < 3; m++) {
        for (int l = 0; l < 3; l++) {
            const LevelMet* L = &MM[m].lk[l];
            printf("%s,%.4f,lock,%d,%d,%.2f,%.2f,%.2f,%.6f,%d,%.3f,%.3f,%.4f\n",
                   MODEL_NAME[m], MM[m].eps, l + 1, L->nv,
                   L->pct[0], L->pct[1], L->pct[2], L->dmax, L->dev_pass,
                   L->uvmax, L->uvrms, L->anchor);
        }
        for (int l = 0; l < 3; l++) {
            const LevelMet* L = &MM[m].fr[l];
            printf("%s,%.4f,free,%d,%d,%.2f,%.2f,%.2f,%.6f,%d,%.3f,%.3f,%.4f\n",
                   MODEL_NAME[m], MM[m].eps, l + 1, L->nv,
                   L->pct[0], L->pct[1], L->pct[2], L->dmax, L->dev_pass,
                   L->uvmax, L->uvrms, L->anchor);
        }
    }
}

int main(int argc, char** argv) {
    for (int m = 0; m < 3; m++) run_model(m);
    if (argc > 1 && strcmp(argv[1], "--csv") == 0) {
        print_csv();
    } else {
        print_report();
    }
    return failures == 0 ? 0 : 1;
}
