
/* symmetry: mirror plane from centroid + principal axis (PCA), density-scaled
 * tolerance, symmetric pair / midline-overlap detection, and ordered cleanup.
 *
 * Pipeline order matters:
 *   raw merged mesh -> detect symmetry + overlaps -> weld -> dedup -> re-count
 * Detecting only AFTER cleanup reads back 0 overlaps and hides the seam defect.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV 96
#define MAXF 128

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { V3 v[MAXV]; int nv; Tri f[MAXF]; int nf; } Mesh;

typedef struct { double x, y, z; } D3;
typedef struct { D3 c; D3 n; double d; double eval[3]; D3 axes[3]; } Plane;

/* density statistics of nearest non-coincident neighbour spacing */
typedef struct { double hmed, hmean, hstd, rho, tol_pair, tol_clean; } Density;

typedef struct {
    int pairs;
    int overlap_pairs_before, overlap_faces_before;
    int overlap_faces_after_weld;
    int overlap_pairs_after, overlap_faces_after;
    int dup_tri_pairs_after_weld;
    double max_mirror_dev;     /* tolerance-free worst nearest-mirror residual */
    double max_pair_residual;  /* worst residual among accepted pairs          */
    int on_plane;
    Plane plane;
    Density dens;
    double bbox_cx;
} Result;

/* ---------------- vector helpers ---------------- */

static D3 d3(double x, double y, double z) { D3 r = { x, y, z }; return r; }
static D3 v3d(V3 p) { return d3(p.x, p.y, p.z); }
static D3 dsub(D3 a, D3 b) { return d3(a.x - b.x, a.y - b.y, a.z - b.z); }
static double ddot(D3 a, D3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static double dnorm(D3 a) { return sqrt(ddot(a, a)); }
static D3 dscale(D3 a, double s) { return d3(a.x*s, a.y*s, a.z*s); }

static double dist3(V3 a, V3 b) {
    double dx = (double)a.x - b.x, dy = (double)a.y - b.y, dz = (double)a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

static int cmp_double(const void* A, const void* B) {
    double a = *(const double*)A, b = *(const double*)B;
    return (a > b) - (a < b);
}

/* ---------------- mesh construction ---------------- */

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

static int add_v(Mesh* m, float x, float y, float z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv] = (V3){ x, y, z };
    return m->nv++;
}

static int add_f(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a < 0 || b < 0 || c < 0 || a >= m->nv || b >= m->nv || c >= m->nv) return -1;
    if (a == b || b == c || a == c) return -1;
    m->f[m->nf] = (Tri){ a, b, c };
    return m->nf++;
}

/* Build a wide "wings" plate from a left half and an independently built right
 * half merged together. x is the span (principal) axis, x=0 the midline.
 * Each side brings its OWN coincident midline vertices and its OWN copy of a
 * midline cap quad, so the raw merge has double vertices + overlapping faces.
 * pairL/pairR (may be NULL) receive the 12 off-plane mirror-pair indices. */
#define NXLEV 3
#define NYROW 4
#define NPAIRS (NXLEV * NYROW)

static int build_wings(Mesh* m, int pairL[NPAIRS], int pairR[NPAIRS]) {
    mesh_init(m);
    const float XMAG[NXLEV] = { 0.9f, 1.8f, 2.7f };
    const float YROW[NYROW] = { 0.0f, 0.75f, 1.5f, 2.25f };
    int base[2][NXLEV + 1][NYROW];
    for (int s = 0; s < 2; s++)
        for (int j = 0; j <= NXLEV; j++)
            for (int r = 0; r < NYROW; r++) base[s][j][r] = -1;

    for (int s = 0; s < 2; s++) {
        float sign = (s == 0) ? -1.0f : 1.0f;   /* 0 = left, 1 = right */
        for (int j = 0; j <= NXLEV; j++) {
            float x = (j == 0) ? 0.0f : sign * XMAG[j - 1];
            for (int r = 0; r < NYROW; r++)
                base[s][j][r] = add_v(m, x, YROW[r], 0.0f);
        }
    }

    if (pairL && pairR) {
        int k = 0;
        for (int j = 1; j <= NXLEV; j++)
            for (int r = 0; r < NYROW; r++) {
                pairL[k] = base[0][j][r];
                pairR[k] = base[1][j][r];
                k++;
            }
    }

    /* skin quads -> 2 triangles each, per side */
    for (int s = 0; s < 2; s++)
        for (int j = 0; j < NXLEV; j++)
            for (int r = 0; r < NYROW - 1; r++) {
                int a = base[s][j][r], b = base[s][j + 1][r];
                int c = base[s][j + 1][r + 1], d = base[s][j][r + 1];
                add_f(m, a, b, c);
                add_f(m, a, c, d);
            }

    /* midline cap quad in the y-z plane, built once per side at x=0.
     * The two copies are geometrically identical -> overlapping midline faces. */
    const float cyz[4][2] = { {0.25f, 0.35f}, {0.25f, -0.35f},
                             {1.85f, -0.30f}, {1.85f, 0.40f} };
    int cap[2][4];
    for (int s = 0; s < 2; s++)
        for (int q = 0; q < 4; q++)
            cap[s][q] = add_v(m, 0.0f, cyz[q][0], cyz[q][1]);
    for (int s = 0; s < 2; s++) {
        add_f(m, cap[s][0], cap[s][1], cap[s][2]);
        add_f(m, cap[s][0], cap[s][2], cap[s][3]);
    }
    return m->nv;
}

/* add a one-sided sparse "spur" at +sx on row y=0.75 with one fin triangle */
static int add_spur(Mesh* m, float sx) {
    int sp = add_v(m, sx, 0.75f, 0.0f);
    /* right outer column (j=3 -> index column NXLEV), rows 1 and 2:
     * vertex layout makes direct lookup awkward, so find by coordinate. */
    int r1 = -1, r2 = -1;
    for (int i = 0; i < m->nv; i++) {
        if (fabsf(m->v[i].x - 2.7f) < 1e-6f && fabsf(m->v[i].y - 0.75f) < 1e-6f) r1 = i;
        if (fabsf(m->v[i].x - 2.7f) < 1e-6f && fabsf(m->v[i].y - 1.50f) < 1e-6f) r2 = i;
    }
    if (sp >= 0 && r1 >= 0 && r2 >= 0) add_f(m, sp, r1, r2);
    return sp;
}

/* ---------------- PCA / mirror plane ---------------- */

static D3 centroid(const Mesh* m) {
    D3 c = d3(0, 0, 0);
    for (int i = 0; i < m->nv; i++) { c.x += m->v[i].x; c.y += m->v[i].y; c.z += m->v[i].z; }
    if (m->nv) c = dscale(c, 1.0 / m->nv);
    return c;
}

/* Jacobi eigen-decomposition of a symmetric 3x3 matrix.
 * On return eval[k] / axes[k] are sorted by descending eigenvalue. */
static void eigen3(double a[3][3], double eval[3], D3 axes[3]) {
    double v[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
    for (int iter = 0; iter < 60; iter++) {
        int p = 0, q = 1;
        double mx = fabs(a[0][1]);
        if (fabs(a[0][2]) > mx) { mx = fabs(a[0][2]); p = 0; q = 2; }
        if (fabs(a[1][2]) > mx) { mx = fabs(a[1][2]); p = 1; q = 2; }
        if (mx < 1e-14) break;
        int r = 3 - p - q;
        double app = a[p][p], aqq = a[q][q], apq = a[p][q];
        double phi = 0.5 * atan2(2.0 * apq, aqq - app);
        double ct = cos(phi), st = sin(phi);
        a[p][p] = ct*ct*app - 2*st*ct*apq + st*st*aqq;
        a[q][q] = st*st*app + 2*st*ct*apq + ct*ct*aqq;
        a[p][q] = a[q][p] = 0.0;
        double arp = a[r][p], arq = a[r][q];
        a[r][p] = a[p][r] = ct*arp - st*arq;
        a[r][q] = a[q][r] = st*arp + ct*arq;
        for (int k = 0; k < 3; k++) {
            double vkp = v[k][p], vkq = v[k][q];
            v[k][p] = ct*vkp - st*vkq;
            v[k][q] = st*vkp + ct*vkq;
        }
    }
    int ord[3] = { 0, 1, 2 };
    double ev[3] = { a[0][0], a[1][1], a[2][2] };
    for (int i = 0; i < 3; i++)
        for (int j = i + 1; j < 3; j++)
            if (ev[ord[j]] > ev[ord[i]]) { int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
    for (int k = 0; k < 3; k++) {
        int col = ord[k];
        eval[k] = ev[col];
        axes[k] = d3(v[0][col], v[1][col], v[2][col]);
        double l = dnorm(axes[k]);
        if (l > 0) axes[k] = dscale(axes[k], 1.0 / l);
    }
}

/* Mirror plane as a function of the vertex centroid and the principal-axis
 * direction: { x : n . x = d }, n = major eigenvector of the covariance,
 * passing through the centroid, hence d = n . c. */
static Plane mirror_plane_pca(const Mesh* m) {
    Plane P; memset(&P, 0, sizeof(P));
    P.c = centroid(m);
    double cov[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
    for (int i = 0; i < m->nv; i++) {
        D3 d = dsub(v3d(m->v[i]), P.c);
        double x[3] = { d.x, d.y, d.z };
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) cov[a][b] += x[a] * x[b];
    }
    if (m->nv > 1)
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) cov[a][b] /= (m->nv - 1);
    eigen3(cov, P.eval, P.axes);
    P.n = P.axes[0];
    P.d = ddot(P.n, P.c);
    return P;
}

static double plane_sd(const Plane* P, V3 p) {
    D3 q = v3d(p);
    return ddot(P->n, q) - P->d;
}

static D3 reflect_point(const Plane* P, V3 p) {
    D3 q = v3d(p);
    double s = ddot(P->n, q) - P->d;
    return dsub(q, dscale(P->n, 2.0 * s));
}

/* ---------------- density-derived tolerance ---------------- */

static double bbox_diag(const Mesh* m) {
    V3 lo = { 1e30f,1e30f,1e30f }, hi = { -1e30f,-1e30f,-1e30f };
    for (int i = 0; i < m->nv; i++) {
        lo.x = fminf(lo.x, m->v[i].x); lo.y = fminf(lo.y, m->v[i].y); lo.z = fminf(lo.z, m->v[i].z);
        hi.x = fmaxf(hi.x, m->v[i].x); hi.y = fmaxf(hi.y, m->v[i].y); hi.z = fmaxf(hi.z, m->v[i].z);
    }
    double dx = hi.x - lo.x, dy = hi.y - lo.y, dz = hi.z - lo.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* Per-vertex nearest NON-coincident neighbour spacing (coincident seam twins
 * are skipped via a bbox-relative floor). Median is robust to the one-sided
 * spur; tolerance is a fixed fraction of that density spacing. */
static Density density_stats(const Mesh* m) {
    Density D; memset(&D, 0, sizeof(D));
    if (m->nv < 2) return D;
    double coin = 1e-6 * fmax(1.0, bbox_diag(m));
    double hs[MAXV]; int nh = 0;
    double sum = 0;
    for (int i = 0; i < m->nv; i++) {
        double best = 1e300;
        for (int j = 0; j < m->nv; j++) {
            if (i == j) continue;
            double dd = dist3(m->v[i], m->v[j]);
            if (dd <= coin) continue;          /* coincident duplicate twin */
            if (dd < best) best = dd;
        }
        if (best < 1e299) { hs[nh++] = best; sum += best; }
    }
    if (nh == 0) return D;
    qsort(hs, nh, sizeof(double), cmp_double);
    D.hmed = (nh % 2) ? hs[nh / 2] : 0.5 * (hs[nh / 2 - 1] + hs[nh / 2]);
    D.hmean = sum / nh;
    double vs = 0;
    for (int i = 0; i < nh; i++) vs += (hs[i] - D.hmean) * (hs[i] - D.hmean);
    D.hstd = sqrt(vs / nh);
    D.rho = 1.0 / (D.hmed * D.hmed);           /* samples per unit area ~ 1/h^2 */
    D.tol_pair = 0.50 * D.hmed;
    D.tol_clean = 0.25 * D.hmed;
    return D;
}

/* ---------------- symmetric pair detection ---------------- */

typedef struct { int li, ri; double res; } Edge;

static int cmp_edge(const void* A, const void* B) {
    double a = ((const Edge*)A)->res, b = ((const Edge*)B)->res;
    return (a > b) - (a < b);
}

/* One-to-one greedy closest assignment among candidate mirror edges inside
 * tol. Returns pair count; *maxres gets worst accepted residual. */
static int detect_pairs(const Mesh* m, const Plane* P, double tol, double* maxres) {
    *maxres = 0.0;
    int neg[MAXV], pos[MAXV], nn = 0, np = 0;
    for (int i = 0; i < m->nv; i++) {
        double s = plane_sd(P, m->v[i]);
        if (s < -tol) neg[nn++] = i;
        else if (s > tol) pos[np++] = i;
    }
    Edge E[MAXV * MAXV]; int ne = 0;
    for (int a = 0; a < nn; a++)
        for (int b = 0; b < np; b++) {
            int li = neg[a], ri = pos[b];
            D3 rp = reflect_point(P, m->v[li]);
            double res = dnorm(dsub(v3d(m->v[ri]), rp));
            if (res <= tol) E[ne++] = (Edge){ li, ri, res };
        }
    qsort(E, ne, sizeof(Edge), cmp_edge);
    int usedL[MAXV] = { 0 }, usedR[MAXV] = { 0 };
    int pairs = 0; double mx = 0;
    for (int k = 0; k < ne; k++) {
        if (usedL[E[k].li] || usedR[E[k].ri]) continue;
        usedL[E[k].li] = usedR[E[k].ri] = 1;
        pairs++;
        if (E[k].res > mx) mx = E[k].res;
    }
    *maxres = mx;
    return pairs;
}

/* tolerance-free worst mirror deviation: for every vertex, distance from its
 * reflected position to the nearest OTHER vertex. Captures unmatched gross
 * asymmetry (spur / lifted verts) instead of hiding it past the pair tol. */
static double max_mirror_deviation(const Mesh* m, const Plane* P) {
    double worst = 0.0;
    for (int i = 0; i < m->nv; i++) {
        D3 rp = reflect_point(P, m->v[i]);
        double best = 1e300;
        for (int j = 0; j < m->nv; j++) {
            if (j == i) continue;
            double dd = dnorm(dsub(v3d(m->v[j]), rp));
            if (dd < best) best = dd;
        }
        if (best > worst) worst = best;
    }
    return worst;
}

static int count_on_plane(const Mesh* m, const Plane* P, double band) {
    int k = 0;
    for (int i = 0; i < m->nv; i++)
        if (fabs(plane_sd(P, m->v[i])) <= band) k++;
    return k;
}

/* ---------------- midline overlap-face detection ---------------- */

static int corner_set_match(const Mesh* m, const Tri* a, const Tri* b, double eps) {
    int ia[3] = { a->a, a->b, a->c };
    int ib[3] = { b->a, b->b, b->c };
    int perm[6][3] = { {0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0} };
    for (int p = 0; p < 6; p++) {
        int ok = 1;
        for (int k = 0; k < 3; k++)
            if (dist3(m->v[ia[k]], m->v[ib[perm[p][k]]]) > eps) { ok = 0; break; }
        if (ok) return 1;
    }
    return 0;
}

/* Count coincident face pairs whose centroid sits on the mirror plane.
 * fpairs (may be NULL) receives pair count; returns involved face count. */
static int detect_midline_overlaps(const Mesh* m, const Plane* P, double eps,
                                   int* fpairs) {
    int paircnt = 0;
    int involved[MAXF]; memset(involved, 0, sizeof(involved));
    for (int i = 0; i < m->nf; i++) {
        for (int j = i + 1; j < m->nf; j++) {
            if (!corner_set_match(m, &m->f[i], &m->f[j], eps)) continue;
            V3 ce = { (m->v[m->f[i].a].x + m->v[m->f[i].b].x + m->v[m->f[i].c].x) / 3.0f,
                      (m->v[m->f[i].a].y + m->v[m->f[i].b].y + m->v[m->f[i].c].y) / 3.0f,
                      (m->v[m->f[i].a].z + m->v[m->f[i].b].z + m->v[m->f[i].c].z) / 3.0f };
            if (fabs(plane_sd(P, ce)) > eps * 2.0 + 1e-9) continue;  /* not midline */
            paircnt++;
            involved[i] = involved[j] = 1;
        }
    }
    if (fpairs) *fpairs = paircnt;
    int faces = 0;
    for (int i = 0; i < m->nf; i++) faces += involved[i];
    return faces;
}

/* count faces sharing the same canonical sorted index triple */
static int count_duplicate_tri_pairs(const Mesh* m) {
    int paircnt = 0;
    for (int i = 0; i < m->nf; i++) {
        int ai[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
        /* sort three ints directly into canonical order */
        for (int x = 0; x < 3; x++) for (int y = x + 1; y < 3; y++)
            if (ai[y] < ai[x]) { int t = ai[x]; ai[x] = ai[y]; ai[y] = t; }
        for (int j = i + 1; j < m->nf; j++) {
            int bi[3] = { m->f[j].a, m->f[j].b, m->f[j].c };
            for (int x = 0; x < 3; x++) for (int y = x + 1; y < 3; y++)
                if (bi[y] < bi[x]) { int t = bi[x]; bi[x] = bi[y]; bi[y] = t; }
            if (ai[0] == bi[0] && ai[1] == bi[1] && ai[2] == bi[2]) paircnt++;
        }
    }
    return paircnt;
}

/* ---------------- ordered cleanup ---------------- */

/* port of welded.c: merge vertices closer than eps, remap, drop degenerates */
static void mesh_weld_eps(Mesh* m, double eps) {
    int map[MAXV];
    V3 uniq[MAXV]; int nu = 0;
    for (int i = 0; i < m->nv; i++) {
        int found = -1;
        for (int j = 0; j < nu; j++)
            if (dist3(m->v[i], uniq[j]) <= eps) { found = j; break; }
        if (found < 0) { uniq[nu] = m->v[i]; found = nu++; }
        map[i] = found;
    }
    int nf = 0;
    for (int i = 0; i < m->nf; i++) {
        int a = map[m->f[i].a], b = map[m->f[i].b], c = map[m->f[i].c];
        if (a == b || b == c || a == c) continue;
        m->f[nf++] = (Tri){ a, b, c };
    }
    m->nf = nf;
    memcpy(m->v, uniq, sizeof(V3) * nu);
    m->nv = nu;
}

/* remove faces with identical canonical index triples; returns removed count */
static int mesh_dedup_faces(Mesh* m) {
    int nf = 0, removed = 0;
    for (int i = 0; i < m->nf; i++) {
        int ai[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
        for (int x = 0; x < 3; x++) for (int y = x + 1; y < 3; y++)
            if (ai[y] < ai[x]) { int t = ai[x]; ai[x] = ai[y]; ai[y] = t; }
        int dup = 0;
        for (int k = 0; k < nf; k++) {
            int bi[3] = { m->f[k].a, m->f[k].b, m->f[k].c };
            for (int x = 0; x < 3; x++) for (int y = x + 1; y < 3; y++)
                if (bi[y] < bi[x]) { int t = bi[x]; bi[x] = bi[y]; bi[y] = t; }
            if (ai[0] == bi[0] && ai[1] == bi[1] && ai[2] == bi[2]) { dup = 1; break; }
        }
        if (dup) { removed++; continue; }
        m->f[nf++] = m->f[i];
    }
    m->nf = nf;
    return removed;
}

static float bbox_cx(const Mesh* m) {
    float lo = 1e30f, hi = -1e30f;
    for (int i = 0; i < m->nv; i++) { lo = fminf(lo, m->v[i].x); hi = fmaxf(hi, m->v[i].x); }
    return 0.5f * (lo + hi);
}

/* ---------------- run one model through the ordered pipeline ---------------- */

static Result analyze(const Mesh* raw) {
    Result R; memset(&R, 0, sizeof(R));
    R.plane = mirror_plane_pca(raw);
    R.dens = density_stats(raw);
    R.bbox_cx = bbox_cx(raw);

    /* 1) DETECT on the raw merged mesh, before any cleanup */
    R.pairs = detect_pairs(raw, &R.plane, R.dens.tol_pair, &R.max_pair_residual);
    R.on_plane = count_on_plane(raw, &R.plane, R.dens.tol_pair);
    R.overlap_faces_before =
        detect_midline_overlaps(raw, &R.plane, R.dens.tol_clean, &R.overlap_pairs_before);
    R.max_mirror_dev = max_mirror_deviation(raw, &R.plane);

    /* 2) weld only: double verts collapse, overlap faces become duplicate tris
     *    but are NOT removed (same behaviour as the welded module). */
    Mesh w = *raw;
    mesh_weld_eps(&w, R.dens.tol_clean);
    R.dup_tri_pairs_after_weld = count_duplicate_tri_pairs(&w);
    R.overlap_faces_after_weld =
        detect_midline_overlaps(&w, &R.plane, R.dens.tol_clean, NULL);

    /* 3) full cleanup: weld + exact-face dedup, then re-count for verification */
    Mesh c = w;
    mesh_dedup_faces(&c);
    R.overlap_faces_after =
        detect_midline_overlaps(&c, &R.plane, R.dens.tol_clean, &R.overlap_pairs_after);
    return R;
}

/* ---------------- self checks ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("check %-28s %s\n", name, ok ? "ok" : "FAIL");
    if (!ok) failures++;
}

static int near(double a, double b, double e) { return fabs(a - b) <= e; }

int main(void) {
    /* ---- unit / invariant checks ---- */
    Mesh t; mesh_init(&t);
    /* axis-aligned anisotropic cross: x is the exact major axis */
    add_v(&t, -2,0,0); add_v(&t, 2,0,0);
    add_v(&t, 0,0.5f,0); add_v(&t, 0,-0.5f,0);
    Plane pt = mirror_plane_pca(&t);
    check("plane_through_centroid", near(pt.d, ddot(pt.n, pt.c), 1e-9));
    check("major_axis_is_x", fabs(fabs(pt.n.x) - 1.0) < 1e-6);

    /* rotated anisotropic cloud: principal axis must follow the distribution */
    Mesh rot; mesh_init(&rot);
    double th = 30.0 * 3.14159265358979 / 180.0, ux = cos(th), uy = sin(th);
    for (int k = -5; k <= 5; k++) {
        double tt = k * 1.0;
        add_v(&rot, (float)(tt*ux), (float)(tt*uy), 0.0f);
        add_v(&rot, (float)(tt*ux + 0.15), (float)(tt*uy - 0.15*ux), 0.2f);
    }
    Plane pr = mirror_plane_pca(&rot);
    double align = fabs(pr.n.x*ux + pr.n.y*uy);
    check("pca_tracks_rotated_axis", align > 0.999);

    /* density spacing scales ~linearly with pitch */
    Mesh g1; mesh_init(&g1);
    for (int x = 0; x < 4; x++) for (int y = 0; y < 4; y++) add_v(&g1, x, y, 0);
    Mesh g2; mesh_init(&g2);
    for (int x = 0; x < 4; x++) for (int y = 0; y < 4; y++) add_v(&g2, 2*x, 2*y, 0);
    Density d1 = density_stats(&g1), d2 = density_stats(&g2);
    check("density_spacing_scales", near(d2.hmed / d1.hmed, 2.0, 0.15));
    check("tol_positive_density_defined", d1.tol_pair > 0 && d1.rho > 0);

    /* sparse one-sided outlier yanks bbox centre but barely moves centroid */
    Mesh ob; mesh_init(&ob);
    for (int k = 0; k < 8; k++) {
        double a = 0.5 + 0.3 * k;
        add_v(&ob, (float)-a, (float)k, 0); add_v(&ob, (float)a, (float)k, 0);
    }
    add_v(&ob, 4.0f, 1.0f, 0.0f);
    check("bbox_center_biased", fabs(bbox_cx(&ob) - (4.0 - 2.6) / 2.0) < 1e-5);
    {
        double ccx = centroid(&ob).x, bx = bbox_cx(&ob);
        check("centroid_far_less_biased", fabs(ccx) < 0.3 && fabs(ccx) < 0.5 * fabs(bx));
    }

    /* ---- three models of graded symmetry ---- */
    Mesh A, B, C;
    int pL[NPAIRS], pR[NPAIRS];

    build_wings(&A, pL, pR);                 /* strict mirror, seam artifacts */

    build_wings(&B, pL, pR);                 /* near-symmetric */
    B.v[pR[2]].z = 0.12f;
    B.v[pR[6]].z = -0.15f;
    B.v[pR[10]].z = 0.18f;
    B.v[pR[5]].z = 0.70f;                    /* one gross lift -> pair lost */
    add_spur(&B, 3.4f);                      /* sparse bbox-biasing spur */

    build_wings(&C, pL, pR);                 /* weakly symmetric */
    C.v[pR[4]].z = 0.90f;                    /* balanced +/- lifts on x=1.8 */
    C.v[pR[5]].z = -0.90f;
    C.v[pR[6]].z = 1.00f;
    C.v[pR[7]].z = -1.00f;                   /* these 4 pairs lost, no PCA tilt */
    C.v[pR[10]].z = 0.90f;                   /* 5th pair lost */
    C.v[pR[1]].z = 0.10f;
    C.v[pR[8]].z = -0.12f;
    add_spur(&C, 4.4f);                      /* large one-sided spur */

    Result rA = analyze(&A), rB = analyze(&B), rC = analyze(&C);

    check("A_all_12_pairs", rA.pairs == NPAIRS);
    check("A_zero_dev", rA.max_mirror_dev < 1e-5);
    check("pairs_monotone", rA.pairs > rB.pairs && rB.pairs > rC.pairs);
    check("dev_monotone", rA.max_mirror_dev < rB.max_mirror_dev &&
                          rB.max_mirror_dev < rC.max_mirror_dev);
    check("overlap_before_present", rA.overlap_faces_before == 4 &&
                                    rB.overlap_faces_before == 4 &&
                                    rC.overlap_faces_before == 4);
    check("weld_alone_keeps_overlap", rA.overlap_faces_after_weld == 4 &&
                                      rA.dup_tri_pairs_after_weld == 2);
    check("dedup_clears_overlap", rA.overlap_faces_after == 0 &&
                                  rB.overlap_faces_after == 0 &&
                                  rC.overlap_faces_after == 0);

    /* ---- detailed per-model dump ---- */
    const char* nm[3] = { "A_strict", "B_near", "C_weak" };
    Result rr[3] = { rA, rB, rC };
    for (int k = 0; k < 3; k++) {
        const Result* R = &rr[k];
        printf("\n=== model %s (nv-plane/density) ===\n", nm[k]);
        printf("centroid=(%.5f,%.5f,%.5f)\n", R->plane.c.x, R->plane.c.y, R->plane.c.z);
        printf("pca_normal=(%.6f,%.6f,%.6f) eig=(%.4f,%.4f,%.6f)\n",
               R->plane.n.x, R->plane.n.y, R->plane.n.z,
               R->plane.eval[0], R->plane.eval[1], R->plane.eval[2]);
        printf("plane_d=%.6f bbox_center_x=%.4f\n", R->plane.d, R->bbox_cx);
        printf("hmed=%.5f hmean=%.5f hstd=%.5f rho=%.4f tol_pair=%.5f tol_clean=%.5f\n",
               R->dens.hmed, R->dens.hmean, R->dens.hstd, R->dens.rho,
               R->dens.tol_pair, R->dens.tol_clean);
        printf("on_plane_verts=%d\n", R->on_plane);
        printf("sym_pairs=%d\n", R->pairs);
        printf("mid_overlap_pairs_before=%d mid_overlap_faces_before=%d\n",
               R->overlap_pairs_before, R->overlap_faces_before);
        printf("max_mirror_dev=%.6f max_pair_residual=%.6f\n",
               R->max_mirror_dev, R->max_pair_residual);
        printf("overlap_faces_after_weld=%d dup_tri_pairs_after_weld=%d\n",
               R->overlap_faces_after_weld, R->dup_tri_pairs_after_weld);
        printf("overlap_pairs_after_clean=%d overlap_faces_after_clean=%d\n",
               R->overlap_pairs_after, R->overlap_faces_after);
    }

    /* ---- the requested three groups / nine numbers ---- */
    printf("\n## NINE: model | sym_pairs | mid_overlap_faces(before clean) | max_mirror_dev\n");
    printf("| model | sym_pairs | mid_overlap_faces | max_mirror_dev |\n");
    printf("|---|---:|---:|---:|\n");
    printf("| A strict | %d | %d | %.6f |\n", rA.pairs, rA.overlap_faces_before, rA.max_mirror_dev);
    printf("| B near   | %d | %d | %.6f |\n", rB.pairs, rB.overlap_faces_before, rB.max_mirror_dev);
    printf("| C weak   | %d | %d | %.6f |\n", rC.pairs, rC.overlap_faces_before, rC.max_mirror_dev);

    /* ---- before/after overlap juxtaposition (ordering evidence) ---- */
    printf("\n## ORDER: midline overlap faces before vs after cleanup\n");
    printf("| model | detect BEFORE clean | after weld-only | after weld+dedup | detect AFTER clean would read |\n");
    printf("|---|---:|---:|---:|---:|\n");
    for (int k = 0; k < 3; k++)
        printf("| %s | %d (%d pairs) | %d | %d | %d |\n", nm[k],
               rr[k].overlap_faces_before, rr[k].overlap_pairs_before,
               rr[k].overlap_faces_after_weld, rr[k].overlap_faces_after,
               rr[k].overlap_faces_after);
    printf("note: weld-only leaves the overlaps as %d.. duplicate tri pairs; "
           "full dedup clears them. Detecting after full clean reports 0 and "
           "hides the raw seam.\n", rA.dup_tri_pairs_after_weld);

    /* ---- bbox-centre vs centroid+PCA plane offset from true x=0 ---- */
    printf("\n## PLANE: naive bbox centre vs centroid+PCA offset (true plane x=0)\n");
    printf("| model | bbox_center_x (off) | pca plane d=n.c (off) |\n|---|---:|---:|\n");
    printf("| A strict | %.4f | %.6f |\n", rA.bbox_cx, rA.plane.d);
    printf("| B near   | %.4f | %.6f |\n", rB.bbox_cx, rB.plane.d);
    printf("| C weak   | %.4f | %.6f |\n", rC.bbox_cx, rC.plane.d);

    /* ---- machine readable ---- */
    printf("\n===JSON===\n");
    printf("{\n");
    printf("  \"tolerance\": {\"statistic\": \"median nearest non-coincident spacing\", "
           "\"pair_factor\": 0.5, \"clean_factor\": 0.25},\n");
    printf("  \"models\": [\n");
    for (int k = 0; k < 3; k++) {
        const Result* R = &rr[k];
        printf("    {\"name\": \"%s\", \"sym_pairs\": %d, "
               "\"mid_overlap_faces_before\": %d, \"max_mirror_dev\": %.6f, "
               "\"overlap_faces_after_weld\": %d, \"overlap_faces_after_clean\": %d, "
               "\"hmed\": %.5f, \"tol_pair\": %.5f, \"bbox_cx\": %.4f, \"pca_d\": %.6f}%s\n",
               nm[k], R->pairs, R->overlap_faces_before, R->max_mirror_dev,
               R->overlap_faces_after_weld, R->overlap_faces_after,
               R->dens.hmed, R->dens.tol_pair, R->bbox_cx, R->plane.d,
               k < 2 ? "," : "");
    }
    printf("  ]\n}\n");

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
