
/* symmetrize: align the mirror plane to the vertex-distribution symmetry axis
 * (centroid + covariance PCA, best-scoring eigen-axis), then average the two
 * sides. Positions are averaged; normals are averaged and re-normalized.
 * Only vertices inside a narrow band around the symmetry plane are touched, so
 * far / genuinely asymmetric parts are preserved.
 *
 * Three procedurally generated character-mask models of decreasing symmetry
 * are run once each. Every number below is computed by this program.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXV 256
#define MAXF 512

typedef struct { double x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

static V3 v3(double x, double y, double z) { V3 r = { x, y, z }; return r; }
static V3 add(V3 a, V3 b) { return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static V3 sub(V3 a, V3 b) { return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static V3 mul(V3 a, double s) { return v3(a.x*s, a.y*s, a.z*s); }
static double dot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static V3 cross(V3 a, V3 b) {
    return v3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}
static double len(V3 a) { return sqrt(dot(a, a)); }
static V3 normz(V3 a) {
    double l = len(a);
    if (l < 1e-30) return v3(0, 0, 0);
    return mul(a, 1.0 / l);
}
static double clampd(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

/* A mirror plane is  n . x = d  (n unit length). */
typedef struct { V3 n; double d; V3 c; } Plane;

static V3 reflect_p(V3 p, Plane pl) {
    double s = dot(pl.n, p) - pl.d;
    return sub(p, mul(pl.n, 2.0 * s));
}
static V3 reflect_v(V3 v, V3 n) {
    return sub(v, mul(n, 2.0 * dot(n, v)));
}

typedef struct {
    V3 p[MAXV];          /* positions */
    V3 n[MAXV];          /* area-weighted vertex normals (unit) */
    int has_n[MAXV];
    Tri f[MAXF];
    int nv, nf;
    V3 true_n; double true_d;   /* ground-truth generation plane */
} Mesh;

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); m->true_n = v3(1,0,0); }

static int add_v(Mesh* m, double x, double y, double z) {
    if (m->nv >= MAXV) return -1;
    m->p[m->nv] = v3(x, y, z);
    return m->nv++;
}
static int add_tri(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a<0||b<0||c<0||a>=m->nv||b>=m->nv||c>=m->nv) return -1;
    if (a==b||b==c||a==c) return -1;
    m->f[m->nf++] = (Tri){ a, b, c };
    return m->nf - 1;
}

/* area-weighted vertex normals (same convention as the normal/ module) */
static void compute_normals(Mesh* m) {
    for (int i = 0; i < m->nv; i++) { m->n[i] = v3(0,0,0); m->has_n[i] = 0; }
    for (int i = 0; i < m->nf; i++) {
        V3 a = m->p[m->f[i].a], b = m->p[m->f[i].b], c = m->p[m->f[i].c];
        V3 fn = cross(sub(b, a), sub(c, a));   /* length = 2*area */
        int ix[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
        for (int k = 0; k < 3; k++) m->n[ix[k]] = add(m->n[ix[k]], fn);
    }
    for (int i = 0; i < m->nv; i++) {
        if (len(m->n[i]) > 1e-20) { m->n[i] = normz(m->n[i]); m->has_n[i] = 1; }
    }
}

/* ---------------- density statistic: median nearest-neighbour spacing ---- */
static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}
static double nn_spacing(const Mesh* m) {
    double vals[MAXV];
    for (int i = 0; i < m->nv; i++) {
        double best = 1e30;
        for (int j = 0; j < m->nv; j++) {
            if (j == i) continue;
            double dd = len(sub(m->p[i], m->p[j]));
            if (dd < best) best = dd;
        }
        vals[i] = best;
    }
    qsort(vals, (size_t)m->nv, sizeof(double), cmp_double);
    return vals[m->nv / 2];
}

/* ---------------- centroid + covariance PCA / Jacobi --------------------- */
typedef struct {
    V3 c;
    double eval[3];
    V3 evec[3];     /* columns / indexed by axis */
} PCA;

/* Jacobi eigendecomposition of a symmetric 3x3 A. On output A's diagonal holds
 * the eigenvalues and V[component][axis] holds the matching eigenvectors. */
static void jacobi3(double A[3][3], double V[3][3]) {
    V[0][0]=1; V[0][1]=0; V[0][2]=0;
    V[1][0]=0; V[1][1]=1; V[1][2]=0;
    V[2][0]=0; V[2][1]=0; V[2][2]=1;
    for (int sweep = 0; sweep < 80; sweep++) {
        double off = fabs(A[0][1]) + fabs(A[0][2]) + fabs(A[1][2]);
        if (off < 1e-15) break;
        for (int p = 0; p < 2; p++) {
            for (int q = p+1; q < 3; q++) {
                if (fabs(A[p][q]) < 1e-18) continue;
                double tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                double t = (tau >= 0 ? 1.0 : -1.0) /
                           (fabs(tau) + sqrt(1.0 + tau*tau));
                double cc = 1.0 / sqrt(1.0 + t*t);
                double ss = t * cc;
                double app = A[p][p], aqq = A[q][q], apq = A[p][q];
                A[p][p] = cc*cc*app - 2*ss*cc*apq + ss*ss*aqq;
                A[q][q] = ss*ss*app + 2*ss*cc*apq + cc*cc*aqq;
                A[p][q] = A[q][p] = 0.0;
                for (int k = 0; k < 3; k++) {
                    if (k == p || k == q) continue;
                    double akp = A[k][p], akq = A[k][q];
                    A[k][p] = A[p][k] = cc*akp - ss*akq;
                    A[k][q] = A[q][k] = ss*akp + cc*akq;
                }
                for (int k = 0; k < 3; k++) {
                    double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = cc*vkp - ss*vkq;
                    V[k][q] = ss*vkp + cc*vkq;
                }
            }
        }
    }
}

static PCA fit_pca(const Mesh* m) {
    PCA pc;
    memset(&pc, 0, sizeof(pc));
    int n = m->nv > 0 ? m->nv : 1;
    V3 c = v3(0,0,0);
    for (int i = 0; i < m->nv; i++) c = add(c, m->p[i]);
    pc.c = mul(c, 1.0 / n);

    double C[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
    for (int i = 0; i < m->nv; i++) {
        V3 d = sub(m->p[i], pc.c);
        double q[3] = { d.x, d.y, d.z };
        for (int a = 0; a < 3; a++)
            for (int b = a; b < 3; b++)
                C[a][b] += q[a]*q[b];
    }
    for (int a = 0; a < 3; a++)
        for (int b = a; b < 3; b++) { C[a][b] /= n; C[b][a] = C[a][b]; }

    double V[3][3];
    jacobi3(C, V);
    for (int a = 0; a < 3; a++) {
        pc.eval[a] = C[a][a];
        pc.evec[a] = normz(v3(V[0][a], V[1][a], V[2][a]));
    }
    return pc;
}

/* mean nearest-reflection residual for a candidate plane, robust median-ish.
 * Self is excluded: otherwise the thin depth axis wins by matching every
 * sheet vertex to itself at 2*|z-zbar| instead of finding its true mirror. */
static double plane_score(const Mesh* m, Plane pl) {
    double vals[MAXV];
    for (int i = 0; i < m->nv; i++) {
        V3 r = reflect_p(m->p[i], pl);
        double best = 1e30;
        for (int j = 0; j < m->nv; j++) {
            if (j == i) continue;
            double dd = len(sub(r, m->p[j]));
            if (dd < best) best = dd;
        }
        vals[i] = best;
    }
    qsort(vals, (size_t)m->nv, sizeof(double), cmp_double);
    return vals[m->nv / 2];   /* median: one-sided outliers do not dominate */
}

/* Offset refinement with the PCA normal held fixed. For mutual reflected
 * neighbours inside an inlier distance, (p_i+p_j)/2 lies on the symmetry
 * plane, so the median of n.midpoint is a robust d. Unlike the centroid this
 * is not pulled by one-sided mass. Orientation is left to centroid+PCA, which
 * is well-conditioned for the wide masks used here; we deliberately do NOT
 * re-fit the normal by ICP, which can fall into a wrong-axis local minimum. */
static Plane refine_plane(const Mesh* m, Plane init, double h) {
    Plane pl = init;
    double inlier = 1.15 * h;
    for (int pass = 0; pass < 2; pass++) {
        int best[MAXV]; double bd[MAXV];
        for (int i = 0; i < m->nv; i++) {
            V3 r = reflect_p(m->p[i], pl);
            best[i] = -1; bd[i] = 1e30;
            for (int j = 0; j < m->nv; j++) {
                if (j == i) continue;
                double dd = len(sub(r, m->p[j]));
                if (dd < bd[i]) { bd[i] = dd; best[i] = j; }
            }
        }
        double dvals[MAXV];
        int nd = 0;
        for (int i = 0; i < m->nv; i++) {
            int j = best[i];
            if (j < 0 || j <= i || best[j] != i) continue;
            if (bd[i] > inlier || bd[j] > inlier) continue;
            V3 mid = mul(add(m->p[i], m->p[j]), 0.5);
            dvals[nd++] = dot(pl.n, mid);
        }
        if (nd < 3) break;
        qsort(dvals, (size_t)nd, sizeof(double), cmp_double);
        pl.d = dvals[nd / 2];
    }
    return pl;
}

/* Fit centroid+PCA, offset-refine every eigen-axis plane, and keep the one
 * with the best mirror score. Refining d BEFORE scoring matters: one-sided
 * mass biases d=n.centroid, adding a systematic 2*offset to the correct axis,
 * which can otherwise let an unbiased-but-wrong axis win the score. */
static Plane fit_mirror_plane(const Mesh* m, int* chosen_axis, PCA* pc_out) {
    PCA pc = fit_pca(m);
    double h = nn_spacing(m);
    int best = 0;
    double best_score = 1e30;
    Plane best_pl = { pc.evec[0], dot(pc.evec[0], pc.c), pc.c };
    for (int a = 0; a < 3; a++) {
        V3 n = pc.evec[a];
        Plane pl = { n, dot(n, pc.c), pc.c };
        pl = refine_plane(m, pl, h);
        double s = plane_score(m, pl);
        if (s < best_score) { best_score = s; best = a; best_pl = pl; }
    }
    Plane pl = best_pl;
    pl.c = pc.c;
    /* stable sign: largest-magnitude component positive; flip d with n */
    int k = 0;
    if (fabs(pl.n.y) > fabs(pl.n.x) && fabs(pl.n.y) >= fabs(pl.n.z)) k = 1;
    if (fabs(pl.n.z) > fabs(pl.n.x) && fabs(pl.n.z) > fabs(pl.n.y)) k = 2;
    double q[3] = { pl.n.x, pl.n.y, pl.n.z };
    if (q[k] < 0) { pl.n = mul(pl.n, -1.0); pl.d = -pl.d; }
    if (chosen_axis) *chosen_axis = best;
    if (pc_out) *pc_out = pc;
    return pl;
}

/* ---------------- metrics ------------------------------------------------ */
typedef struct { double mean, mx; } ErrStat;

static ErrStat mirror_error(const Mesh* m, Plane pl, const int* ids, int nid) {
    double sum = 0, mx = 0;
    for (int a = 0; a < nid; a++) {
        int i = ids[a];
        V3 r = reflect_p(m->p[i], pl);
        double best = 1e30;
        for (int j = 0; j < m->nv; j++) {
            double dd = len(sub(r, m->p[j]));
            if (dd < best) best = dd;
        }
        if (best > mx) mx = best;
        sum += best;
    }
    ErrStat e = { nid ? sum / nid : 0.0, mx };
    return e;
}

static double angle_deg(V3 a, V3 b) {
    double c = dot(normz(a), normz(b));
    return acos(clampd(c, -1.0, 1.0)) * 180.0 / 3.14159265358979323846;
}

/* ---------------- symmetrize -------------------------------------------- */
typedef struct {
    int roi_n;
    int processed_v;
    double band, tol;
    ErrStat err_pre, err_post;     /* against operating plane, ROI vertices */
    ErrStat true_pre, true_post;   /* against generation plane, ROI vertices */
    double ang_pre_mean, ang_post_mean, ang_post_max;
    double min_n_len_after;
    double raw_avg_n_len;          /* shortest un-renormalized averaged normal */
    int n_pairs;
} SymResult;

/* Symmetrize in place towards plane pl. band/tol are in distance units. */
static SymResult symmetrize(Mesh* m, Plane pl, double band, double tol) {
    SymResult R;
    memset(&R, 0, sizeof(R));
    R.band = band; R.tol = tol; R.min_n_len_after = 1.0; R.raw_avg_n_len = 1.0;

    int inband[MAXV] = { 0 };
    int roi_ids[MAXV] = { 0 }; int nroi = 0;
    for (int i = 0; i < m->nv; i++) {
        double s = dot(pl.n, m->p[i]) - pl.d;
        if (fabs(s) <= band) { inband[i] = 1; roi_ids[nroi++] = i; }
    }
    R.roi_n = nroi;
    R.err_pre = mirror_error(m, pl, roi_ids, nroi);
    {
        Plane tp = { m->true_n, m->true_d, v3(0,0,0) };
        R.true_pre = mirror_error(m, tp, roi_ids, nroi);
    }

    /* best reflected counterpart for every vertex, from ORIGINAL positions */
    int best[MAXV]; double bd[MAXV];
    for (int i = 0; i < m->nv; i++) {
        V3 r = reflect_p(m->p[i], pl);
        best[i] = -1; bd[i] = 1e30;
        for (int j = 0; j < m->nv; j++) {
            double dd = len(sub(r, m->p[j]));
            if (dd < bd[i]) { bd[i] = dd; best[i] = j; }
        }
    }

    int done[MAXV] = { 0 };
    double apre = 0, apost = 0, apostmax = 0; int nang = 0;
    for (int i = 0; i < m->nv; i++) {
        int j = best[i];
        if (j <= i || done[i] || done[j]) continue;
        if (best[j] != i) continue;                 /* mutual pair only */
        if (bd[i] > tol || bd[j] > tol) continue;
        if (!inband[i] || !inband[j]) continue;     /* only near-plane band */

        /* pre normal disagreement: n_i vs reflected n_j */
        double pre = 0;
        V3 rnj = m->has_n[j] ? reflect_v(m->n[j], pl.n) : v3(0,0,0);
        if (m->has_n[i] && m->has_n[j]) pre = angle_deg(m->n[i], rnj);

        /* average positions across the plane */
        V3 mi = mul(add(m->p[i], reflect_p(m->p[j], pl)), 0.5);
        V3 mj = reflect_p(mi, pl);

        /* average normals in i-frame, then re-normalize (both sides) */
        V3 ni2 = m->n[i], nj2 = rnj;
        V3 avg = mul(add(ni2, nj2), 0.5);
        double rawlen = len(avg);
        if (rawlen < R.raw_avg_n_len) R.raw_avg_n_len = rawlen;
        V3 ani = normz(avg);
        V3 anj = reflect_v(ani, pl.n);

        m->p[i] = mi; m->p[j] = mj;
        if (m->has_n[i] && m->has_n[j]) {
            m->n[i] = ani; m->n[j] = normz(anj);
        }
        done[i] = done[j] = 1;
        R.processed_v += 2; R.n_pairs++;

        if (m->has_n[i] && m->has_n[j]) {
            apre += pre;
            V3 chk_i = m->n[i], chk_j = reflect_v(m->n[j], pl.n);
            double post = angle_deg(chk_i, chk_j);
            apost += post; if (post > apostmax) apostmax = post; nang++;
            double li = len(m->n[i]), lj = len(m->n[j]);
            if (li < R.min_n_len_after) R.min_n_len_after = li;
            if (lj < R.min_n_len_after) R.min_n_len_after = lj;
        }
    }
    if (nang) {
        R.ang_pre_mean = apre / nang;
        R.ang_post_mean = apost / nang;
        R.ang_post_max = apostmax;
    }
    R.err_post = mirror_error(m, pl, roi_ids, nroi);
    {
        Plane tp = { m->true_n, m->true_d, v3(0,0,0) };
        R.true_post = mirror_error(m, tp, roi_ids, nroi);
    }
    return R;
}

/* ---------------- procedural character-mask models ----------------------- */
typedef struct {
    const char* key;
    const char* label;
    int severity;         /* scales the far-feature asymmetry */
    double tx, tilt;      /* rigid placement versus the canonical plane */
} ModelCfg;

/* Build an oval relief mask. Left/right grid columns are exact mirrors until
 * we perturb right-side central vertices and add asymmetric far features. */
static Mesh build_model(ModelCfg cfg, Plane* true_plane_world) {
    Mesh m; mesh_init(&m);
    enum { NX = 8, NY = 9 };
    double dx = 0.92, dy = 0.45;
    int idx[NX][NY];

    for (int ix = 0; ix < NX; ix++) {
        double x = (ix - (NX-1)/2.0) * dx;
        for (int iy = 0; iy < NY; iy++) {
            double y = iy * dy;
            /* face oval: narrow the jaw (low y) and forehead (high y) */
            double yn = (y - 4*dy) / (4*dy);          /* -1..1 */
            double profile = 0.72 + 0.28 * (1.0 - yn*yn);
            double xx = x * profile;
            /* relief: nose bump off the vertical centre so the model is not
             * also symmetric top/bottom (which would confuse PCA axis pick) */
            double nose = 0.42 * exp(-(xx*xx)/0.5) *
                          exp(-((y-1.55)*(y-1.55))/0.55);
            double brow = 0.10 * exp(-((y-2.5)*(y-2.5))/0.5);
            double z = nose + brow;
            idx[ix][iy] = add_v(&m, xx, y, z);
        }
    }
    for (int ix = 0; ix < NX-1; ix++)
        for (int iy = 0; iy < NY-1; iy++) {
            int a = idx[ix][iy], b = idx[ix+1][iy];
            int c = idx[ix+1][iy+1], d = idx[ix][iy+1];
            add_tri(&m, a, b, c); add_tri(&m, a, c, d);
        }

    /* Explicit detail warps on the innermost RIGHT column ix=4 (x~+0.28),
     * which is inside the symmetry band. Offsets alternate sign so the
     * covariance cross-term stays small and PCA keeps the true axis; large
     * warps exceed the pairing tolerance and are deliberately left unpaired. */
    if (strcmp(cfg.key, "mid") == 0) {
        m.p[idx[4][2]].y += 0.18;
        m.p[idx[4][3]].z += 0.16;
        m.p[idx[5][5]].y -= 0.17;
        m.p[idx[4][6]].z += 0.20;
        m.p[idx[4][5]].y += 0.46;   /* just over tol: this pair is skipped */
    } else if (strcmp(cfg.key, "low") == 0) {
        /* depth (z) warps give mirror/normal error but do not move vertices
         * out of the near-plane band (the plane normal has ~no z) and do not
         * torque the xy PCA axis. 3 oversized pairs are refused and keep a
         * large residual; 4 moderate pairs are still matched and averaged. */
        m.p[idx[4][1]].z += 0.24;
        m.p[idx[4][2]].z += 0.68;   /* refused */
        m.p[idx[4][3]].z += 0.30;
        m.p[idx[4][4]].z -= 0.74;   /* refused */
        m.p[idx[4][5]].z += 0.27;
        m.p[idx[4][6]].z -= 0.66;   /* refused */
        m.p[idx[4][7]].z += 0.33;
    }

    /* Far feature triangles: exactly 6 vertices / 2 triangles for every real
     * model so denominators match. strict = mirrored horns, mid = one drooped
     * horn, low = two +x spikes placed above/below centroid so their xy
     * covariance torque cancels (axis stays recoverable). All far outside band.
     * The "clean" PCA fixture adds none. */
    if (strcmp(cfg.key, "clean") != 0) {
        double sev = 0.10 * cfg.severity;
        int h1 = m.nv, h2;
        if (strcmp(cfg.key, "low") == 0) {
            add_v(&m, 2.20, 2.45, 0.05);
            add_v(&m, 2.55 + sev, 2.55, 0.02);
            add_v(&m, 2.24, 2.86, 0.10);
            add_tri(&m, h1, h1+1, h1+2);
            h2 = m.nv;
            add_v(&m, 2.24, 1.30, 0.05);
            add_v(&m, 2.60 + sev, 1.20, 0.02);
            add_v(&m, 2.26, 0.92, 0.10);
            add_tri(&m, h2, h2+1, h2+2);
        } else {
            /* shared +x upper horn */
            add_v(&m, 2.45, 3.30, 0.05);
            add_v(&m, 2.78, 3.36, 0.02);
            add_v(&m, 2.52, 3.74, 0.12);
            add_tri(&m, h1, h1+1, h1+2);
            h2 = m.nv;
            if (strcmp(cfg.key, "strict") == 0) {
                add_v(&m, -2.45, 3.30, 0.05);
                add_v(&m, -2.78, 3.36, 0.02);
                add_v(&m, -2.52, 3.74, 0.12);
            } else { /* mid: near-mirror horn, slightly drooped */
                double dr = 0.16 + 0.05 * cfg.severity;
                add_v(&m, -2.45, 3.30 - dr, 0.05);
                add_v(&m, -2.78, 3.36 - dr, 0.02);
                add_v(&m, -2.52, 3.74 - dr, 0.12);
            }
            add_tri(&m, h2, h2+2, h2+1);
        }
    }

    /* rigid placement: translate x and tilt around z (lean in xy plane) */
    double al = cfg.tilt;
    double ca = cos(al), sa = sin(al);
    for (int i = 0; i < m.nv; i++) {
        double x = m.p[i].x, y = m.p[i].y;
        m.p[i].x = ca*x - sa*y + cfg.tx;
        m.p[i].y = sa*x + ca*y;
    }
    /* true plane after transform: normal = Rz(tilt)*(1,0,0) = (cos,sin,0),
     * through transformed origin point (tx,0,0). */
    V3 tn = v3(ca, sa, 0);
    V3 p0 = v3(cfg.tx, 0, 0);
    m.true_n = tn; m.true_d = dot(tn, p0);
    if (true_plane_world) *true_plane_world = (Plane){ tn, m.true_d, p0 };

    compute_normals(&m);
    return m;
}

/* ---------------- top-level experiment ----------------------------------- */
typedef struct {
    ModelCfg cfg;
    Plane fit; int axis; PCA pca;
    double h;
    SymResult sr;
    int total_v;
    double pct;
} Row;

static Row run_aligned(ModelCfg cfg) {
    Row r; memset(&r, 0, sizeof(r)); r.cfg = cfg;
    Plane tp;
    Mesh m = build_model(cfg, &tp);
    r.total_v = m.nv;
    r.h = nn_spacing(&m);
    r.fit = fit_mirror_plane(&m, &r.axis, &r.pca);
    double B = 1.55 * r.h;
    double tol = 0.85 * r.h;
    r.sr = symmetrize(&m, r.fit, B, tol);
    r.pct = 100.0 * r.sr.processed_v / r.total_v;
    return r;
}

/* forced misaligned plane (e.g. naive world x=0), operate on fresh mesh */
static Row run_forced_plane(ModelCfg cfg, Plane bad) {
    Row r; memset(&r, 0, sizeof(r)); r.cfg = cfg;
    Plane tp;
    Mesh m = build_model(cfg, &tp);
    r.total_v = m.nv;
    r.h = nn_spacing(&m);
    r.fit = bad;
    double B = 1.55 * r.h;
    double tol = 0.85 * r.h;
    r.sr = symmetrize(&m, bad, B, tol);
    r.pct = 100.0 * r.sr.processed_v / r.total_v;
    return r;
}

/* The exact diagnostic predicate from the brief:
 * least-symmetry row has the highest processed fraction AND its objective
 * (true-plane) mean symmetry error did not improve. */
static int misaligned_signature(const Row rows[3]) {
    const Row* lo = &rows[2];
    int highest = 1;
    for (int k = 0; k < 3; k++)
        if (rows[k].pct > lo->pct + 1e-9) highest = 0;
    double pre = lo->sr.true_pre.mean, post = lo->sr.true_post.mean;
    int no_gain = post >= pre - 1e-6;     /* no strict improvement */
    return highest && no_gain;
}

/* ============================ self tests ================================= */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int approx(double a, double b, double e) { return fabs(a-b) <= e; }

static void self_tests(void) {
    /* reflection math */
    Plane xp = { v3(1,0,0), 0, v3(0,0,0) };
    V3 q = reflect_p(v3(1,2,3), xp);
    check("reflect_point", approx(q.x,-1,1e-12)&&approx(q.y,2,1e-12)&&approx(q.z,3,1e-12));
    V3 rv = reflect_v(v3(1,0,0), v3(1,0,0));
    check("reflect_vec", approx(rv.x,-1,1e-12));

    /* PCA recovers a known tilted, offset plane */
    {
        ModelCfg c = { "clean","clean tilted fixture", 0, 0.7, 0.18 };
        Plane tp; Mesh m = build_model(c, &tp);
        int ax; PCA pc; Plane fit = fit_mirror_plane(&m, &ax, &pc);
        double ang = angle_deg(fit.n, tp.n);
        if (ang > 90) ang = 180 - ang;      /* normals sign-equivalent */
        check("pca_tilt_axis", ang < 2.0);
        check("pca_offset_plane", fabs(fit.d - tp.d) < 0.25);
        check("centroid_reported", len(sub(fit.c, pc.c)) < 1e-9);
    }

    /* exact symmetric model -> near zero error and full inner pairing */
    {
        ModelCfg c = { "strict","A", 1, 0.0, 0.0 };
        Row r = run_aligned(c);
        check("strict_pre_small", r.sr.err_pre.mean < 1e-6);
        check("strict_post_zero", r.sr.err_post.mx < 1e-6);
        check("strict_normals_agree", r.sr.ang_post_max < 1e-4);
        check("normals_unit_after", fabs(r.sr.min_n_len_after-1.0) < 1e-9);
        check("pairs_were_done", r.sr.n_pairs > 0);
    }

    /* un-renormalized averaged normal shrinks when normals disagree: this is
     * the collapse that re-normalization prevents */
    {
        V3 a = v3(0,0,1);
        V3 b = normz(v3(0.35,0,0.94));      /* ~20 deg apart */
        V3 avg = mul(add(a,b),0.5);
        check("raw_normal_shrinks", len(avg) < 0.999);
        check("renorm_restores_unit", fabs(len(normz(avg))-1.0) < 1e-12);
    }

    /* averaging normals is required: without it paired normals stay apart */
    {
        ModelCfg c = { "mid","B", 1, 0.0, 0.0 };
        Plane tp; Mesh m = build_model(c,&tp);
        Plane fit = fit_mirror_plane(&m,NULL,NULL);
        SymResult full = symmetrize(&m, fit, 1.55*nn_spacing(&m), 0.85*nn_spacing(&m));
        check("normal_blend_zeroes_angle", full.ang_post_mean < 1e-6);
        /* the small jitter must actually show up as pre-normal disagreement */
        check("normal_pre_disagreement_seen", full.ang_pre_mean > 1.0);
    }

    /* far asymmetric vertices are never moved */
    {
        ModelCfg c = { "low","C", 3, 0.0, 0.0 };
        Plane tp; Mesh m = build_model(c,&tp);
        Mesh before = m;
        Plane fit = fit_mirror_plane(&m,NULL,NULL);
        double h = nn_spacing(&m);
        /* identify far vertices by |signed| > band and require unchanged */
        double B = 1.55*h;
        int moved_far = 0;
        /* symmetrize mutates; snapshot then compare */
        V3 oldp[MAXV]; memcpy(oldp, m.p, sizeof(oldp));
        symmetrize(&m, fit, B, 0.85*h);
        for (int i=0;i<before.nv;i++){
            double s = dot(fit.n, oldp[i]) - fit.d;
            if (fabs(s) > B + 1e-9 && len(sub(m.p[i], oldp[i])) > 1e-12) moved_far++;
        }
        check("far_vertices_untouched", moved_far == 0);
    }

    /* fitted plane must recover a tilted/offset weak-symmetry model, and the
     * operation must improve its objective symmetry, with coverage ordered
     * strict > mid > low */
    {
        ModelCfg sc = { "strict","A", 1, 0.00, 0.00 };
        ModelCfg mc = { "mid",   "B", 1, 0.10, 0.05 };
        ModelCfg lc = { "low",   "C", 3, 0.34, 0.16 };
        Row rs = run_aligned(sc), rm = run_aligned(mc), rl = run_aligned(lc);
        Plane tp; Mesh lm = build_model(lc, &tp); (void)lm;
        double ang = angle_deg(rl.fit.n, tp.n);
        if (ang > 90) ang = 180 - ang;
        check("low_plane_axis_recovered", ang < 2.0);
        check("low_plane_offset_recovered", fabs(rl.fit.d - tp.d) < 0.06);
        check("low_objective_improves",
              rl.sr.true_post.mean < rl.sr.true_pre.mean - 0.02);
        check("coverage_monotone", rs.pct > rm.pct && rm.pct > rl.pct);
        check("error_monotone",
              rs.sr.err_post.mx < rm.sr.err_post.mx &&
              rm.sr.err_post.mx < rl.sr.err_post.mx);
    }

    /* diagnostic predicate: aligned healthy rows do NOT flag; a constructed
     * misaligned triplet that matches the brief DOES flag. */
    {
        ModelCfg cfgs[3] = {
            { "strict","A", 1, 0, 0 },
            { "mid",   "B", 1, 0, 0 },
            { "low",   "C", 3, 0, 0 },
        };
        Row rows[3];
        for (int k=0;k<3;k++) rows[k]=run_aligned(cfgs[k]);
        check("aligned_no_false_alarm", misaligned_signature(rows)==0);

        Row bad[3];
        for (int k=0;k<3;k++) { bad[k]=rows[k]; bad[k].pct = 20.0 + 5.0*k; }
        bad[2].pct = 60.0;                       /* low now highest */
        bad[2].sr.true_pre.mean = 0.40;
        bad[2].sr.true_post.mean = 0.43;        /* and got worse */
        check("misaligned_predicate_fires", misaligned_signature(bad)==1);

        Row bad2[3];
        for (int k=0;k<3;k++) bad2[k]=rows[k];
        check("healthy_real_rows_not_flagged", misaligned_signature(bad2)==0);
    }
}

/* ============================ reporting ================================== */
static const ModelCfg CFGS[3] = {
    { "strict", "A strict (near-perfect mirror)", 1, 0.00, 0.00 },
    { "mid",    "B approximate (small jitter)",    1, 0.10, 0.05 },
    { "low",    "C weak (large one-sided warp)",   3, 0.34, 0.16 },
};

static void print_plane(const char* indent, Plane pl) {
    printf("%scentroid=(%.6f, %.6f, %.6f)  plane: %.6f*x %+.6f*y %+.6f*z = %.6f"
           "   (d=n.c=%.6f)\n",
           indent, pl.c.x, pl.c.y, pl.c.z,
           pl.n.x, pl.n.y, pl.n.z, pl.d, dot(pl.n, pl.c));
}

static void print_human(Row rows[3], Row* wrong) {
    printf("\n================ SYMMETRIZE: 3 MODELS ================\n");
    for (int k=0;k<3;k++) {
        Row* r=&rows[k];
        printf("\n[%s] %s  nv=%d h=%.4f band=%.3f tol=%.3f chosen_eigen_axis=%d\n",
               r->cfg.key, r->cfg.label, r->total_v, r->h,
               r->sr.band, r->sr.tol, r->axis);
        printf("  eigenvalues=%.5g %.5g %.5g\n",
               r->pca.eval[0], r->pca.eval[1], r->pca.eval[2]);
        print_plane("  fitted ", r->fit);
        {
            Plane tp; Mesh m=build_model(r->cfg,&tp); (void)m;
            printf("  true    plane: %.6f*x %+.6f*y %+.6f*z = %.6f\n",
                   tp.n.x,tp.n.y,tp.n.z,tp.d);
        }
        printf("  ROI vertices=%d  processed=%d (%.2f%% of all %d) pairs=%d\n",
               r->sr.roi_n, r->sr.processed_v, r->pct, r->total_v, r->sr.n_pairs);
        printf("  mirror err vs fit-plane  mean %.6f -> %.6f   max %.6f -> %.6f\n",
               r->sr.err_pre.mean, r->sr.err_post.mean,
               r->sr.err_pre.mx, r->sr.err_post.mx);
        printf("  objective err vs true    mean %.6f -> %.6f   max %.6f -> %.6f\n",
               r->sr.true_pre.mean, r->sr.true_post.mean,
               r->sr.true_pre.mx, r->sr.true_post.mx);
        printf("  normal angle deg(pre mean %.4f) -> post mean %.6f max %.6f | "
               "min|n| after %.9f | raw avg min|n| %.6f\n",
               r->sr.ang_pre_mean, r->sr.ang_post_mean, r->sr.ang_post_max,
               r->sr.min_n_len_after, r->sr.raw_avg_n_len);
    }

    printf("\n---------------- three groups, nine numbers ----------------\n");
    printf("canonical nine = 3 models x [err_mean_post, err_max_post, normal_angle_post_deg]\n");
    printf("model, err_mean_post, err_max_post, normal_angle_post_deg\n");
    for (int k=0;k<3;k++)
        printf("NINE_ROW %s %.6f %.6f %.6f\n", rows[k].cfg.key,
               rows[k].sr.err_post.mean, rows[k].sr.err_post.mx,
               rows[k].sr.ang_post_mean);
    printf("NINE_VECTOR");
    for (int k=0;k<3;k++)
        printf(" %.6f %.6f %.6f", rows[k].sr.err_post.mean,
               rows[k].sr.err_post.mx, rows[k].sr.ang_post_mean);
    printf("\n");

    printf("\nfull record: symmetry error pre/post, normal angle pre/post, coverage\n");
    printf("model errMeanPre errMeanPost errMaxPre errMaxPost "
           "angPreDeg angPostDeg processedPct\n");
    for (int k=0;k<3;k++)
        printf("DETAIL_ROW %s %.6f %.6f %.6f %.6f %.4f %.6f %.2f\n",
               rows[k].cfg.key,
               rows[k].sr.err_pre.mean, rows[k].sr.err_post.mean,
               rows[k].sr.err_pre.mx, rows[k].sr.err_post.mx,
               rows[k].sr.ang_pre_mean, rows[k].sr.ang_post_mean,
               rows[k].pct);

    printf("\n---------------- alignment diagnostic ----------------\n");
    printf("aligned-plane misalignment_signature=%d (expect 0)\n",
           misaligned_signature(rows));
    if (wrong) {
        printf("\nCONTROL: same low-symmetry model forced onto a WRONG "
               "horizontal axis plane n=(0,1,0) through centroid y.\n");
        print_plane("  fitted(correct) ", rows[2].fit);
        print_plane("  forced(wrong)   ", wrong[2].fit);
        printf("  processed %.2f%% (correct) -> %.2f%% (wrong)\n",
               rows[2].pct, wrong[2].pct);
        printf("  objective true-plane err mean %.6f -> %.6f (wrong op)\n",
               wrong[2].sr.true_pre.mean, wrong[2].sr.true_post.mean);
        printf("  wrong-plane self residual     %.6f -> %.6f\n",
               wrong[2].sr.err_pre.mean, wrong[2].sr.err_post.mean);
        printf("  CONTROL low processed higher than aligned: %d ; "
               "objective error improved: %d\n",
               wrong[2].pct > rows[2].pct + 1e-9,
               wrong[2].sr.true_post.mean < wrong[2].sr.true_pre.mean - 1e-6);

        /* 3-row forced-wrong table for the literal predicate */
        int hi=1; for (int k=0;k<3;k++) if (wrong[k].pct>wrong[2].pct+1e-9) hi=0;
        int nogain = wrong[2].sr.true_post.mean >= wrong[2].sr.true_pre.mean-1e-6;
        printf("\nforced-wrong 3-row processed%%: %.2f %.2f %.2f\n",
               wrong[0].pct, wrong[1].pct, wrong[2].pct);
        printf("forced-wrong low objective mean err pre/post: %.6f / %.6f\n",
               wrong[2].sr.true_pre.mean, wrong[2].sr.true_post.mean);
        printf("forced-wrong signature components: low_is_max_pct=%d no_gain=%d "
               "=> misaligned_signature=%d\n", hi, nogain,
               misaligned_signature(wrong));
        printf("FITTED_PLANE_VS_CENTROID low: c=(%.6f,%.6f,%.6f) "
               "A=%.8f B=%.8f C=%.8f D=%.8f\n",
               rows[2].fit.c.x,rows[2].fit.c.y,rows[2].fit.c.z,
               rows[2].fit.n.x,rows[2].fit.n.y,rows[2].fit.n.z,rows[2].fit.d);
    }
}

static void print_json(Row rows[3], Row wrong[3]) {
    printf("{\n  \"experiment\": \"symmetrize: PCA mirror-plane alignment + "
           "position/normal averaging, near-plane band only\",\n");
    printf("  \"reproduce\": \"make -C symmetrize test\",\n");
    printf("  \"nine_order\": \"row-major, each model = [err_mean_post, "
           "err_max_post, normal_angle_post_deg]\",\n");
    printf("  \"nine\": [");
    int first=1;
    for (int k=0;k<3;k++) {
        double v[3]={rows[k].sr.err_post.mean,rows[k].sr.err_post.mx,
                     rows[k].sr.ang_post_mean};
        for (int t=0;t<3;t++){ if(!first) printf(","); first=0; printf(" %.6f",v[t]); }
    }
    printf(" ],\n  \"models\": [\n");
    for (int k=0;k<3;k++) {
        Row* r=&rows[k];
        printf("    {\"key\":\"%s\",\"nv\":%d,\"h\":%.6f,\"band\":%.6f,\"tol\":%.6f,"
               "\"chosen_axis\":%d,\"centroid\":[%.6f,%.6f,%.6f],"
               "\"plane_n\":[%.8f,%.8f,%.8f],\"plane_d\":%.8f,"
               "\"plane_d_through_centroid\":%.8f,"
               "\"eigenvalues\":[%.6f,%.6f,%.6f],"
               "\"roi_n\":%d,\"processed_v\":%d,\"processed_pct\":%.4f,"
               "\"err_mean_pre\":%.6f,\"err_max_pre\":%.6f,"
               "\"err_mean_post\":%.6f,\"err_max_post\":%.6f,"
               "\"true_mean_pre\":%.6f,\"true_max_pre\":%.6f,"
               "\"true_mean_post\":%.6f,\"true_max_post\":%.6f,"
               "\"normal_angle_pre_mean\":%.6f,\"normal_angle_post_mean\":%.8f,"
               "\"normal_angle_post_max\":%.8f,\"min_normal_len\":%.9f}%s\n",
               r->cfg.key,r->total_v,r->h,r->sr.band,r->sr.tol,r->axis,
               r->fit.c.x,r->fit.c.y,r->fit.c.z,
               r->fit.n.x,r->fit.n.y,r->fit.n.z,r->fit.d,
               dot(r->fit.n, r->fit.c),
               r->pca.eval[0],r->pca.eval[1],r->pca.eval[2],
               r->sr.roi_n,r->sr.processed_v,r->pct,
               r->sr.err_pre.mean,r->sr.err_pre.mx,
               r->sr.err_post.mean,r->sr.err_post.mx,
               r->sr.true_pre.mean,r->sr.true_pre.mx,
               r->sr.true_post.mean,r->sr.true_post.mx,
               r->sr.ang_pre_mean,r->sr.ang_post_mean,r->sr.ang_post_max,
               r->sr.min_n_len_after, k<2?",":"");
    }
    printf("  ],\n");
    printf("  \"aligned_misaligned_signature\": %d,\n", misaligned_signature(rows));
    printf("  \"control_wrong_axis\": [\n");
    for (int k=0;k<3;k++) {
        Row* r=&wrong[k];
        printf("    {\"key\":\"%s\",\"processed_pct\":%.4f,"
               "\"true_mean_pre\":%.6f,\"true_mean_post\":%.6f,"
               "\"self_mean_pre\":%.6f,\"self_mean_post\":%.6f}%s\n",
               r->cfg.key,r->pct,r->sr.true_pre.mean,r->sr.true_post.mean,
               r->sr.err_pre.mean,r->sr.err_post.mean, k<2?",":"");
    }
    printf("  ],\n  \"control_misaligned_signature\": %d\n}\n",
           misaligned_signature(wrong));
}

int main(int argc, char** argv) {
    int json = (argc > 1 && strcmp(argv[1], "--json") == 0);

    Row rows[3], wrong[3];
    for (int k=0;k<3;k++) rows[k] = run_aligned(CFGS[k]);

    /* forced wrong horizontal plane through each model centroid */
    for (int k=0;k<3;k++) {
        Plane bad = { v3(0,1,0), rows[k].fit.c.y, rows[k].fit.c };
        wrong[k] = run_forced_plane(CFGS[k], bad);
    }

    if (json) {
        print_json(rows, wrong);
        return 0;
    }

    self_tests();
    print_human(rows, wrong);
    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
