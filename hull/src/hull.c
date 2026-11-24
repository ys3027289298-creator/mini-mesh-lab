
/* hull: 3D convex hull of a triangle mesh's vertex set, no external deps.
 *
 * The hull of a mesh equals the hull of its vertex set, so face
 * connectivity is not needed on input. Two independent implementations
 * share the same face-pool machinery:
 *   - quickhull:   seed tetra from extreme points, then repeatedly expand
 *                  at the furthest outside point of some face
 *   - incremental: seed tetra from the first non-degenerate 4 points in
 *                  input order, then insert points one by one
 * Both emit the hull as (vertex indices into input, triangles, volume),
 * so vertex sets and volumes can be cross-checked between the two.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

typedef struct { double x, y, z; } V3;

static V3 v3sub(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static V3 v3cross(V3 a, V3 b) {
    V3 r = { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
    return r;
}
static double v3dot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static double v3len(V3 a) { return sqrt(v3dot(a, a)); }

/* ---------- face pool ---------- */

typedef struct {
    int v[3];     /* input point indices, CCW seen from outside */
    V3 n;         /* outward unit normal */
    double d;     /* plane: n.x = d */
    int alive;
    int* pts;     /* quickhull only: outside points assigned to this face */
    int npts, cap;
} Face;

typedef struct { Face* f; int n, cap; } FPool;

static void fp_init(FPool* p) { p->f = NULL; p->n = 0; p->cap = 0; }

static void fp_free(FPool* p) {
    for (int i = 0; i < p->n; i++) free(p->f[i].pts);
    free(p->f);
    fp_init(p);
}

static double face_dist(const Face* f, V3 p) { return v3dot(f->n, p) - f->d; }

/* add face through (a,b,c); flip winding so `inside` lies behind the plane */
static int fp_add(FPool* p, const V3* v, int a, int b, int c, V3 inside) {
    if (p->n == p->cap) {
        int nc = p->cap ? p->cap * 2 : 16;
        Face* nf = realloc(p->f, sizeof(Face) * nc);
        if (!nf) return -1;
        p->f = nf; p->cap = nc;
    }
    V3 e1 = v3sub(v[b], v[a]);
    V3 e2 = v3sub(v[c], v[a]);
    V3 nrm = v3cross(e1, e2);
    double L = v3len(nrm);
    if (L == 0.0) return -1;
    nrm.x /= L; nrm.y /= L; nrm.z /= L;
    if (v3dot(nrm, v3sub(inside, v[a])) > 0.0) {
        int t = b; b = c; c = t;
        nrm.x = -nrm.x; nrm.y = -nrm.y; nrm.z = -nrm.z;
    }
    Face* f = &p->f[p->n];
    f->v[0] = a; f->v[1] = b; f->v[2] = c;
    f->n = nrm;
    f->d = v3dot(nrm, v[a]);
    f->alive = 1;
    f->pts = NULL; f->npts = 0; f->cap = 0;
    return p->n++;
}

static void face_add_pt(Face* f, int idx) {
    if (f->npts == f->cap) {
        int nc = f->cap ? f->cap * 2 : 8;
        f->pts = realloc(f->pts, sizeof(int) * nc);
        f->cap = nc;
    }
    f->pts[f->npts++] = idx;
}

static int face_has_vert(const Face* f, int g) {
    return f->v[0] == g || f->v[1] == g || f->v[2] == g;
}

/* mark faces visible from point pi; emit horizon edges (pairs of input
 * indices). vis must hold fp->n chars, hor must hold 3*fp->n pairs.
 * returns the horizon edge count (0 = point not outside any face). */
static int collect_horizon(const FPool* fp, const V3* v, int pi, double eps,
                           char* vis, int (*hor)[2]) {
    for (int i = 0; i < fp->n; i++)
        vis[i] = (fp->f[i].alive && face_dist(&fp->f[i], v[pi]) > eps) ? 1 : 0;
    int nh = 0;
    for (int i = 0; i < fp->n; i++) {
        if (!vis[i]) continue;
        const Face* f = &fp->f[i];
        for (int e = 0; e < 3; e++) {
            int a = f->v[e], b = f->v[(e + 1) % 3];
            int adj_visible = 0;
            for (int g = 0; g < fp->n; g++) {
                if (g == i || !fp->f[g].alive) continue;
                if (face_has_vert(&fp->f[g], a) && face_has_vert(&fp->f[g], b) && vis[g]) {
                    adj_visible = 1;
                    break;
                }
            }
            if (!adj_visible) { hor[nh][0] = a; hor[nh][1] = b; nh++; }
        }
    }
    return nh;
}

/* ---------- hull result ---------- */

typedef struct {
    int* idx;      /* hull vertex indices into input, first-appearance order */
    int nv;
    int (*tri)[3]; /* faces as indices into idx */
    int nf;
    double volume;
} Hull;

static void hull_free(Hull* h) {
    free(h->idx);
    free(h->tri);
    memset(h, 0, sizeof(*h));
}

static int hull_extract(Hull* h, const FPool* fp, const V3* v) {
    int nf = 0;
    for (int i = 0; i < fp->n; i++) if (fp->f[i].alive) nf++;
    if (nf == 0) return -1;
    h->idx = NULL; h->nv = 0;
    h->tri = malloc(sizeof(int[3]) * nf);
    h->nf = nf;
    int cap = 0, k = 0;
    double vol = 0.0;
    for (int i = 0; i < fp->n; i++) {
        if (!fp->f[i].alive) continue;
        const Face* f = &fp->f[i];
        for (int j = 0; j < 3; j++) {
            int g = f->v[j], pos = -1;
            for (int q = 0; q < h->nv; q++) if (h->idx[q] == g) { pos = q; break; }
            if (pos < 0) {
                if (h->nv == cap) {
                    cap = cap ? cap * 2 : 16;
                    h->idx = realloc(h->idx, sizeof(int) * cap);
                }
                pos = h->nv;
                h->idx[h->nv++] = g;
            }
            h->tri[k][j] = pos;
        }
        vol += v3dot(v[f->v[0]], v3cross(v[f->v[1]], v[f->v[2]]));
        k++;
    }
    h->volume = vol / 6.0;
    return 0;
}

/* ---------- seed tetrahedra ---------- */

static int seed_tetra_extreme(const V3* p, int n, int t[4]) {
    int ex[6] = { 0, 0, 0, 0, 0, 0 };
    for (int i = 1; i < n; i++) {
        if (p[i].x < p[ex[0]].x) ex[0] = i;
        if (p[i].x > p[ex[1]].x) ex[1] = i;
        if (p[i].y < p[ex[2]].y) ex[2] = i;
        if (p[i].y > p[ex[3]].y) ex[3] = i;
        if (p[i].z < p[ex[4]].z) ex[4] = i;
        if (p[i].z > p[ex[5]].z) ex[5] = i;
    }
    double best = -1.0;
    int a = 0, b = 1;
    for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++) {
            V3 d = v3sub(p[ex[i]], p[ex[j]]);
            double d2 = v3dot(d, d);
            if (d2 > best) { best = d2; a = ex[i]; b = ex[j]; }
        }
    if (best == 0.0) return -1; /* all points coincide */
    t[0] = a; t[1] = b;
    V3 ab = v3sub(p[b], p[a]);
    double L = v3len(ab);
    best = -1.0;
    int c = -1;
    for (int i = 0; i < n; i++) {
        double d = v3len(v3cross(v3sub(p[i], p[a]), ab)) / L;
        if (d > best) { best = d; c = i; }
    }
    if (c < 0 || best == 0.0) return -1; /* collinear */
    t[2] = c;
    V3 nrm = v3cross(v3sub(p[t[1]], p[t[0]]), v3sub(p[t[2]], p[t[0]]));
    best = -1.0;
    int dd = -1;
    for (int i = 0; i < n; i++) {
        double d = fabs(v3dot(nrm, v3sub(p[i], p[t[0]])));
        if (d > best) { best = d; dd = i; }
    }
    if (dd < 0 || best == 0.0) return -1; /* coplanar */
    t[3] = dd;
    return 0;
}

static int seed_tetra_seq(const V3* p, int n, int t[4]) {
    t[0] = 0;
    t[1] = t[2] = t[3] = -1;
    for (int i = 1; i < n; i++)
        if (v3len(v3sub(p[i], p[t[0]])) > 0.0) { t[1] = i; break; }
    if (t[1] < 0) return -1;
    V3 ab = v3sub(p[t[1]], p[t[0]]);
    for (int i = t[1] + 1; i < n; i++)
        if (v3len(v3cross(ab, v3sub(p[i], p[t[0]]))) > 0.0) { t[2] = i; break; }
    if (t[2] < 0) return -1;
    V3 nrm = v3cross(ab, v3sub(p[t[2]], p[t[0]]));
    for (int i = t[2] + 1; i < n; i++)
        if (fabs(v3dot(nrm, v3sub(p[i], p[t[0]]))) > 0.0) { t[3] = i; break; }
    if (t[3] < 0) return -1;
    return 0;
}

static V3 tetra_centroid(const V3* p, const int t[4]) {
    V3 c = { 0, 0, 0 };
    for (int i = 0; i < 4; i++) { c.x += p[t[i]].x; c.y += p[t[i]].y; c.z += p[t[i]].z; }
    c.x /= 4.0; c.y /= 4.0; c.z /= 4.0;
    return c;
}

static double auto_eps(const V3* p, int n) {
    V3 lo = p[0], hi = p[0];
    for (int i = 1; i < n; i++) {
        if (p[i].x < lo.x) lo.x = p[i].x;
        if (p[i].y < lo.y) lo.y = p[i].y;
        if (p[i].z < lo.z) lo.z = p[i].z;
        if (p[i].x > hi.x) hi.x = p[i].x;
        if (p[i].y > hi.y) hi.y = p[i].y;
        if (p[i].z > hi.z) hi.z = p[i].z;
    }
    double diag = v3len(v3sub(hi, lo));
    return 1e-9 * (diag > 0.0 ? diag : 1.0);
}

/* ---------- algorithm 1: quickhull ---------- */

static int hull_quickhull(const V3* p, int n, double eps, Hull* out) {
    if (n < 4) return -1;
    if (eps <= 0.0) eps = auto_eps(p, n);
    int t[4];
    if (seed_tetra_extreme(p, n, t) != 0) return -1;
    V3 inside = tetra_centroid(p, t);
    FPool fp;
    fp_init(&fp);
    static const int cf[4][3] = { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3} };
    for (int i = 0; i < 4; i++)
        fp_add(&fp, p, t[cf[i][0]], t[cf[i][1]], t[cf[i][2]], inside);
    for (int i = 0; i < n; i++) {
        if (i == t[0] || i == t[1] || i == t[2] || i == t[3]) continue;
        for (int f = 0; f < fp.n; f++)
            if (face_dist(&fp.f[f], p[i]) > eps) { face_add_pt(&fp.f[f], i); break; }
    }
    int guard = 2 * n + 16;
    while (guard-- > 0) {
        int f = -1;
        for (int i = 0; i < fp.n; i++)
            if (fp.f[i].alive && fp.f[i].npts > 0) { f = i; break; }
        if (f < 0) break;
        int pi = -1;
        double best = eps;
        for (int j = 0; j < fp.f[f].npts; j++) {
            double d = face_dist(&fp.f[f], p[fp.f[f].pts[j]]);
            if (d > best) { best = d; pi = fp.f[f].pts[j]; }
        }
        if (pi < 0) { fp.f[f].npts = 0; continue; }
        char* vis = malloc(fp.n);
        int (*hor)[2] = malloc(sizeof(int[2]) * 3 * (fp.n + 1));
        int nh = collect_horizon(&fp, p, pi, eps, vis, hor);
        if (nh == 0) { /* numeric noise: drop pi from its list */
            for (int j = 0; j < fp.f[f].npts; j++)
                if (fp.f[f].pts[j] == pi) { fp.f[f].pts[j] = fp.f[f].pts[--fp.f[f].npts]; break; }
            free(vis); free(hor);
            continue;
        }
        int cap = 0;
        for (int i = 0; i < fp.n; i++) if (vis[i]) cap += fp.f[i].npts;
        int* orph = malloc(sizeof(int) * (cap > 0 ? cap : 1));
        int no = 0;
        for (int i = 0; i < fp.n; i++)
            if (vis[i])
                for (int j = 0; j < fp.f[i].npts; j++)
                    if (fp.f[i].pts[j] != pi) orph[no++] = fp.f[i].pts[j];
        for (int i = 0; i < fp.n; i++)
            if (vis[i]) { fp.f[i].alive = 0; free(fp.f[i].pts); fp.f[i].pts = NULL; fp.f[i].npts = 0; }
        int first_new = fp.n;
        for (int e = 0; e < nh; e++)
            fp_add(&fp, p, hor[e][0], hor[e][1], pi, inside);
        for (int o = 0; o < no; o++)
            for (int f2 = first_new; f2 < fp.n; f2++)
                if (face_dist(&fp.f[f2], p[orph[o]]) > eps) { face_add_pt(&fp.f[f2], orph[o]); break; }
        free(vis); free(hor); free(orph);
    }
    int rc = (guard <= 0) ? -1 : hull_extract(out, &fp, p);
    fp_free(&fp);
    return rc;
}

/* ---------- algorithm 2: point-by-point incremental ---------- */

static int hull_incremental(const V3* p, int n, double eps, Hull* out) {
    if (n < 4) return -1;
    if (eps <= 0.0) eps = auto_eps(p, n);
    int t[4];
    if (seed_tetra_seq(p, n, t) != 0) return -1;
    V3 inside = tetra_centroid(p, t);
    FPool fp;
    fp_init(&fp);
    static const int cf[4][3] = { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3} };
    for (int i = 0; i < 4; i++)
        fp_add(&fp, p, t[cf[i][0]], t[cf[i][1]], t[cf[i][2]], inside);
    char* seeded = calloc(n, 1);
    for (int i = 0; i < 4; i++) seeded[t[i]] = 1;
    for (int i = 0; i < n; i++) {
        if (seeded[i]) continue;
        char* vis = malloc(fp.n);
        int (*hor)[2] = malloc(sizeof(int[2]) * 3 * (fp.n + 1));
        int nh = collect_horizon(&fp, p, i, eps, vis, hor);
        if (nh > 0) {
            for (int j = 0; j < fp.n; j++) if (vis[j]) fp.f[j].alive = 0;
            for (int e = 0; e < nh; e++)
                fp_add(&fp, p, hor[e][0], hor[e][1], i, inside);
        }
        free(vis); free(hor);
    }
    free(seeded);
    int rc = hull_extract(out, &fp, p);
    fp_free(&fp);
    return rc;
}

/* ---------- validation helpers ---------- */

/* every input point must lie behind (or on) every hull face */
static int hull_contains_all(const Hull* h, const V3* p, int n, double tol) {
    for (int f = 0; f < h->nf; f++) {
        V3 a = p[h->idx[h->tri[f][0]]];
        V3 b = p[h->idx[h->tri[f][1]]];
        V3 c = p[h->idx[h->tri[f][2]]];
        V3 nrm = v3cross(v3sub(b, a), v3sub(c, a));
        double L = v3len(nrm);
        if (L == 0.0) return 0;
        nrm.x /= L; nrm.y /= L; nrm.z /= L;
        double d = v3dot(nrm, a);
        for (int i = 0; i < n; i++)
            if (v3dot(nrm, p[i]) - d > tol) return 0;
    }
    return 1;
}

static int cmp_int(const void* a, const void* b) {
    int x = *(const int*)a, y = *(const int*)b;
    return (x > y) - (x < y);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* cross-check the two algorithms on the same model */
static void compare_hulls(const char* name, const Hull* qh, const Hull* inc) {
    int* sa = malloc(sizeof(int) * qh->nv);
    int* sb = malloc(sizeof(int) * inc->nv);
    memcpy(sa, qh->idx, sizeof(int) * qh->nv);
    memcpy(sb, inc->idx, sizeof(int) * inc->nv);
    qsort(sa, qh->nv, sizeof(int), cmp_int);
    qsort(sb, inc->nv, sizeof(int), cmp_int);
    int set_eq = (qh->nv == inc->nv);
    if (set_eq)
        for (int i = 0; i < qh->nv; i++)
            if (sa[i] != sb[i]) { set_eq = 0; break; }
    int first_diff = -1;
    int m = qh->nv < inc->nv ? qh->nv : inc->nv;
    for (int i = 0; i < m; i++)
        if (qh->idx[i] != inc->idx[i]) { first_diff = i; break; }
    if (first_diff < 0 && qh->nv != inc->nv) first_diff = m;
    double dv = fabs(qh->volume - inc->volume);
    double ref = fabs(qh->volume) > 1.0 ? fabs(qh->volume) : 1.0;
    char buf[64];
    snprintf(buf, sizeof buf, "cmp_%s_set_equal", name); check(buf, set_eq);
    snprintf(buf, sizeof buf, "cmp_%s_vol_match", name); check(buf, dv / ref < 1e-9);
    printf("cmp_%s_order_first_diff=%d\n", name, first_diff);
    printf("cmp_%s_vol_qh=%.15g\n", name, qh->volume);
    printf("cmp_%s_vol_inc=%.15g\n", name, inc->volume);
    printf("cmp_%s_vol_rel_diff=%.3g\n", name, dv / ref);
    printf("cmp_%s_order_qh=", name);
    for (int i = 0; i < qh->nv && i < 8; i++) printf("%d ", qh->idx[i]);
    printf("\n");
    printf("cmp_%s_order_inc=", name);
    for (int i = 0; i < inc->nv && i < 8; i++) printf("%d ", inc->idx[i]);
    printf("\n");
    free(sa); free(sb);
}

/* ---------- test model generators ---------- */

static int gen_sphere(V3* v, int nseg, int nring) {
    int n = 0;
    v[n++] = (V3){ 0, 0, 1 };
    for (int r = 1; r < nring; r++) {
        double phi = PI * r / nring;
        for (int s = 0; s < nseg; s++) {
            double th = 2.0 * PI * s / nseg;
            v[n++] = (V3){ sin(phi) * cos(th), sin(phi) * sin(th), cos(phi) };
        }
    }
    v[n++] = (V3){ 0, 0, -1 };
    return n;
}

/* cube [-1,1]^3 + 4x4 grid on the top face with the 4 inner grid points
 * pushed down to z=0.5 (the dent) + 3 strictly interior points */
static int gen_dented_cube(V3* v) {
    int n = 0;
    for (int x = 0; x < 2; x++)
        for (int y = 0; y < 2; y++)
            for (int z = 0; z < 2; z++)
                v[n++] = (V3){ x ? 1.0 : -1.0, y ? 1.0 : -1.0, z ? 1.0 : -1.0 };
    double g[4] = { -1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0 };
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            if ((i == 0 || i == 3) && (j == 0 || j == 3)) continue; /* corners added */
            double z = (i == 0 || i == 3 || j == 0 || j == 3) ? 1.0 : 0.5;
            v[n++] = (V3){ g[i], g[j], z };
        }
    v[n++] = (V3){ 0, 0, 0 };
    v[n++] = (V3){ 0.2, -0.3, 0.1 };
    v[n++] = (V3){ -0.4, 0.1, -0.2 };
    return n;
}

int main(void) {
    /* 1. single tetrahedron: hull is itself, V = 1/6 */
    V3 tet[4] = { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} };
    Hull qt, it;
    check("tet_qh_ok", hull_quickhull(tet, 4, 0, &qt) == 0);
    check("tet_qh_nv", qt.nv == 4);
    check("tet_qh_nf", qt.nf == 4);
    check("tet_qh_vol", fabs(qt.volume - 1.0 / 6.0) < 1e-12);
    check("tet_qh_contains", hull_contains_all(&qt, tet, 4, 1e-9));
    check("tet_inc_ok", hull_incremental(tet, 4, 0, &it) == 0);
    check("tet_inc_nv", it.nv == 4);
    check("tet_inc_vol", fabs(it.volume - 1.0 / 6.0) < 1e-12);

    /* 2. hollow sphere: every vertex lies on the sphere, hull keeps all */
    static V3 sph[512];
    int ns = gen_sphere(sph, 24, 12);
    Hull qs, is;
    check("sphere_nv", ns == 266);
    check("sphere_qh_ok", hull_quickhull(sph, ns, 0, &qs) == 0);
    check("sphere_inc_ok", hull_incremental(sph, ns, 0, &is) == 0);
    check("sphere_qh_all_kept", qs.nv == ns);
    check("sphere_inc_all_kept", is.nv == ns);
    double vsph = 4.0 / 3.0 * PI;
    check("sphere_qh_vol", qs.volume < vsph && qs.volume > 0.9 * vsph);
    check("sphere_inc_vol", is.volume < vsph && is.volume > 0.9 * vsph);
    check("sphere_qh_contains", hull_contains_all(&qs, sph, ns, 1e-6));
    check("sphere_inc_contains", hull_contains_all(&is, sph, ns, 1e-6));

    /* 3. dented cube: hull collapses to the 8 corners, V = 8 */
    V3 cube[64];
    int nc = gen_dented_cube(cube);
    Hull qc, ic;
    check("cube_input_nv", nc == 23);
    check("cube_qh_ok", hull_quickhull(cube, nc, 0, &qc) == 0);
    check("cube_inc_ok", hull_incremental(cube, nc, 0, &ic) == 0);
    check("cube_qh_nv", qc.nv == 8);
    check("cube_inc_nv", ic.nv == 8);
    check("cube_qh_vol", fabs(qc.volume - 8.0) < 1e-9);
    check("cube_inc_vol", fabs(ic.volume - 8.0) < 1e-9);
    check("cube_qh_contains", hull_contains_all(&qc, cube, nc, 1e-9));
    check("cube_inc_contains", hull_contains_all(&ic, cube, nc, 1e-9));
    check("cube_hull_much_smaller", qc.nv * 2 < nc);

    /* 4. cross-check both algorithms on the same models */
    compare_hulls("cube", &qc, &ic);
    compare_hulls("sphere", &qs, &is);

    hull_free(&qt); hull_free(&it);
    hull_free(&qs); hull_free(&is);
    hull_free(&qc); hull_free(&ic);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
