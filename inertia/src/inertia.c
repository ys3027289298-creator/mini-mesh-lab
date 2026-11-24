/* inertia: volume, centroid, inertia tensor and principal axes of a closed
 * triangle mesh, from a per-tetrahedron (origin, a, b, c) decomposition.
 *
 * The raw moments P = integral(x x^T dV) and m1 = integral(x dV) come out of
 * the tet decomposition referred to the ORIGIN.  Everything below about the
 * centroid is obtained by the parallel-axis shift
 *
 *     S   = P - V c c^T                      (second moment about the centroid)
 *     I_c = tr(S) I - S = I_o - V (|c|^2 I - c c^T)
 *
 * The shape principal axes are the eigenvectors of S, the inertia principal
 * axes those of I_c.  Both are computed here; the "origin" tensor is also
 * reported so the caller can see what the un-shifted decomposition would have
 * produced.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

typedef struct { double x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct {
    V3 *v; int nv, capv;
    Tri *f; int nf, capf;
} Mesh;

typedef struct {
    double V;
    V3 centroid;
    double P[3][3];        /* integral x_i x_j dV about the origin */
    double S[3][3];        /* second moment about the centroid */
    double I_origin[3][3]; /* tr(P) I - P: what you get without the shift   */
    double I[3][3];        /* tr(S) I - S: inertia about the centroid       */
    int closed;
} Moments;

/* ---------- vectors / matrices ---------- */

static V3 v3(double x, double y, double z) { V3 r = {x, y, z}; return r; }
static V3 cross(V3 p, V3 q) {
    return v3(p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x);
}
static double dot(V3 p, V3 q) { return p.x * q.x + p.y * q.y + p.z * q.z; }
static double norm(V3 p) { return sqrt(dot(p, p)); }

static void mat_zero(double A[3][3]) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) A[i][j] = 0.0;
}

static void mat_identity(double A[3][3]) {
    mat_zero(A);
    A[0][0] = A[1][1] = A[2][2] = 1.0;
}

/* inertia tensor from a second-moment matrix about the same point */
static void inertia_from_second_moment(const double S[3][3], double I[3][3]) {
    double tr = S[0][0] + S[1][1] + S[2][2];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            I[i][j] = (i == j ? tr : 0.0) - S[i][j];
}

/* ---------- mesh ---------- */

static void mesh_init(Mesh *m) { memset(m, 0, sizeof(*m)); }
static void mesh_free(Mesh *m) { free(m->v); free(m->f); mesh_init(m); }

static int mesh_add_vertex(Mesh *m, V3 p) {
    if (m->nv == m->capv) {
        int cap = m->capv ? m->capv * 2 : 16;
        V3 *nv = realloc(m->v, sizeof(V3) * (size_t)cap);
        if (!nv) return -1;
        m->v = nv; m->capv = cap;
    }
    m->v[m->nv++] = p;
    return 0;
}

static int mesh_add_tri(Mesh *m, int a, int b, int c) {
    if (m->nf == m->capf) {
        int cap = m->capf ? m->capf * 2 : 16;
        Tri *nf = realloc(m->f, sizeof(Tri) * (size_t)cap);
        if (!nf) return -1;
        m->f = nf; m->capf = cap;
    }
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    m->nf++;
    return 0;
}

/* add an axis-aligned box as 12 outward-facing (CCW) triangles */
static void mesh_add_box(Mesh *m, double x0, double y0, double z0, double x1, double y1, double z1) {
    int base = m->nv;
    mesh_add_vertex(m, v3(x0, y0, z0)); /* 0 */
    mesh_add_vertex(m, v3(x1, y0, z0)); /* 1 */
    mesh_add_vertex(m, v3(x1, y1, z0)); /* 2 */
    mesh_add_vertex(m, v3(x0, y1, z0)); /* 3 */
    mesh_add_vertex(m, v3(x0, y0, z1)); /* 4 */
    mesh_add_vertex(m, v3(x1, y0, z1)); /* 5 */
    mesh_add_vertex(m, v3(x1, y1, z1)); /* 6 */
    mesh_add_vertex(m, v3(x0, y1, z1)); /* 7 */
    static const int t[12][3] = {
        {0, 2, 1}, {0, 3, 2}, /* -z */
        {4, 5, 6}, {4, 6, 7}, /* +z */
        {0, 1, 5}, {0, 5, 4}, /* -y */
        {2, 3, 7}, {2, 7, 6}, /* +y */
        {0, 4, 7}, {0, 7, 3}, /* -x */
        {1, 2, 6}, {1, 6, 5}, /* +x */
    };
    for (int i = 0; i < 12; i++)
        mesh_add_tri(m, base + t[i][0], base + t[i][1], base + t[i][2]);
}

/* ---------- OBJ input ---------- */

static char *trim(char *s) {
    while (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n') s++;
    char *e = s + strlen(s);
    while (e > s && (e[-1] == ' ' || e[-1] == '\t' || e[-1] == '\r' || e[-1] == '\n')) *--e = '\0';
    return s;
}

static int parse_obj(const char *path, Mesh *m) {
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        char *s = trim(line);
        if (s[0] == 'v' && (s[1] == ' ' || s[1] == '\t')) {
            double x, y, z;
            if (sscanf(s + 1, "%lf %lf %lf", &x, &y, &z) == 3)
                if (mesh_add_vertex(m, v3(x, y, z)) != 0) { fclose(fp); return -1; }
        } else if (s[0] == 'f' && (s[1] == ' ' || s[1] == '\t')) {
            int idx[16], n = 0;
            char *p = s + 1;
            while (n < 16) {
                while (*p == ' ' || *p == '\t') p++;
                if (!*p) break;
                char *end;
                long k = strtol(p, &end, 10);
                if (end == p) break;
                if (k < 0) k = m->nv + k + 1;      /* 1-based, negatives from end */
                idx[n++] = (int)k - 1;
                p = end;
                while (*p && *p != ' ' && *p != '\t') p++; /* skip /vt/vn */
            }
            for (int i = 1; i + 1 < n; i++)
                if (mesh_add_tri(m, idx[0], idx[i], idx[i + 1]) != 0) { fclose(fp); return -1; }
        }
    }
    fclose(fp);
    return 0;
}

/* ---------- topology check ---------- */

typedef struct { int u, v, dir; } Edge;

static int edge_cmp(const void *pa, const void *pb) {
    const Edge *a = pa, *b = pb;
    if (a->u != b->u) return a->u - b->u;
    if (a->v != b->v) return a->v - b->v;
    return a->dir - b->dir;
}

/* closed + consistently oriented: every undirected edge is used once in each
 * direction.  The origin-tet formulas are only valid on such a mesh. */
static int mesh_is_closed(const Mesh *m) {
    if (m->nf == 0) return 0;
    Edge *e = malloc(sizeof(Edge) * (size_t)m->nf * 3);
    if (!e) return 0;
    int n = 0;
    for (int i = 0; i < m->nf; i++) {
        int tri[3] = {m->f[i].a, m->f[i].b, m->f[i].c};
        for (int k = 0; k < 3; k++) {
            int a = tri[k], b = tri[(k + 1) % 3];
            e[n].u = a < b ? a : b;
            e[n].v = a < b ? b : a;
            e[n].dir = a < b ? 1 : -1;
            n++;
        }
    }
    qsort(e, (size_t)n, sizeof(Edge), edge_cmp);
    int ok = 1;
    for (int i = 0; i < n;) {
        int j = i, sum = 0;
        while (j < n && e[j].u == e[i].u && e[j].v == e[i].v) { sum += e[j].dir; j++; }
        if (j - i != 2 || sum != 0) { ok = 0; break; }
        i = j;
    }
    free(e);
    return ok;
}

/* ---------- moments ---------- */

static void mesh_moments(const Mesh *m, Moments *mom) {
    double V6 = 0.0;
    V3 m1 = v3(0, 0, 0);
    double P[3][3];
    mat_zero(P);

    for (int i = 0; i < m->nf; i++) {
        V3 a = m->v[m->f[i].a], b = m->v[m->f[i].b], c = m->v[m->f[i].c];
        double det = dot(a, cross(b, c));       /* 6 * signed tet volume */
        V6 += det;
        m1 = v3(m1.x + det * (a.x + b.x + c.x),
                m1.y + det * (a.y + b.y + c.y),
                m1.z + det * (a.z + b.z + c.z));

        /* integral of x x^T over the tet (0,a,b,c):
         *   det/120 * (sum_i p_i p_i^T + (a+b+c)(a+b+c)^T),  det = 6 dV */
        double pv[3][3] = {{a.x, a.y, a.z}, {b.x, b.y, b.z}, {c.x, c.y, c.z}};
        double u[3] = {a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z};
        for (int r = 0; r < 3; r++)
            for (int s = 0; s < 3; s++) {
                double acc = 0.0;
                for (int k = 0; k < 3; k++) acc += pv[k][r] * pv[k][s];
                P[r][s] += det / 120.0 * (acc + u[r] * u[s]);
            }
    }

    mom->V = V6 / 6.0;
    mom->centroid = v3(m1.x / 24.0 / mom->V, m1.y / 24.0 / mom->V, m1.z / 24.0 / mom->V);

    /* parallel-axis shift from the origin to the centroid */
    double c[3] = {mom->centroid.x, mom->centroid.y, mom->centroid.z};
    for (int r = 0; r < 3; r++)
        for (int s = 0; s < 3; s++) {
            mom->P[r][s] = P[r][s];
            mom->S[r][s] = P[r][s] - mom->V * c[r] * c[s];
        }

    inertia_from_second_moment(mom->P, mom->I_origin);
    inertia_from_second_moment(mom->S, mom->I);

    mom->closed = mesh_is_closed(m);
}

/* ---------- symmetric 3x3 eigen-decomposition (cyclic Jacobi) ---------- */

static void jacobi3(const double Ain[3][3], double eval[3], double evec[3][3]) {
    double A[3][3], V[3][3];
    memcpy(A, Ain, sizeof(A));
    mat_identity(V);

    for (int sweep = 0; sweep < 100; sweep++) {
        double off = fabs(A[0][1]) + fabs(A[0][2]) + fabs(A[1][2]);
        double diag = fabs(A[0][0]) + fabs(A[1][1]) + fabs(A[2][2]);
        if (off <= 1e-17 * diag || off == 0.0) break;
        for (int p = 0; p < 2; p++) {
            for (int q = p + 1; q < 3; q++) {
                if (A[p][q] == 0.0) continue;
                double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                double t = (theta >= 0.0 ? 1.0 : -1.0) /
                           (fabs(theta) + sqrt(theta * theta + 1.0));
                double cs = 1.0 / sqrt(t * t + 1.0);
                double sn = t * cs;
                double app = A[p][p], aqq = A[q][q], apq = A[p][q];
                A[p][p] = cs * cs * app - 2.0 * sn * cs * apq + sn * sn * aqq;
                A[q][q] = sn * sn * app + 2.0 * sn * cs * apq + cs * cs * aqq;
                A[p][q] = A[q][p] = 0.0;
                for (int k = 0; k < 3; k++) {
                    if (k == p || k == q) continue;
                    double akp = A[k][p], akq = A[k][q];
                    A[k][p] = A[p][k] = cs * akp - sn * akq;
                    A[k][q] = A[q][k] = sn * akp + cs * akq;
                }
                for (int k = 0; k < 3; k++) {
                    double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = cs * vkp - sn * vkq;
                    V[k][q] = sn * vkp + cs * vkq;
                }
            }
        }
    }
    for (int i = 0; i < 3; i++) {
        eval[i] = A[i][i];
        for (int k = 0; k < 3; k++) evec[i][k] = V[k][i]; /* row i = eigenvector i */
    }
}

/* sort eigenpairs by eigenvalue; ascending for inertia, descending for shape */
static void sort_eigenpairs(double eval[3], double evec[3][3], int ascending) {
    for (int i = 0; i < 2; i++) {
        for (int j = i + 1; j < 3; j++) {
            bool swap = ascending ? (eval[j] < eval[i]) : (eval[j] > eval[i]);
            if (swap) {
                double t = eval[i]; eval[i] = eval[j]; eval[j] = t;
                for (int k = 0; k < 3; k++) {
                    t = evec[i][k]; evec[i][k] = evec[j][k]; evec[j][k] = t;
                }
            }
        }
    }
}

/* angle in degrees between two eigenvectors, sign-insensitive */
static double axis_angle_deg(const double a[3], const double b[3]) {
    const double rad_to_deg = 180.0 / 3.14159265358979323846;
    double d = fabs(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
    if (d > 1.0) d = 1.0;
    return acos(d) * rad_to_deg;
}

/* Assign each shape axis (already in descending-second-moment order) to the
 * inertia axis that is closest to it.  Six permutations, pick the one with the
 * smallest total angle; this stays stable when two eigenvalues are close. */
static void match_axes(const double shape[3][3], const double inertia[3][3],
                       double angles[3], int match[3]) {
    static const int perm[6][3] = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    double best = 1e300;
    for (int p = 0; p < 6; p++) {
        double total = 0.0;
        for (int i = 0; i < 3; i++) total += axis_angle_deg(shape[i], inertia[perm[p][i]]);
        if (total < best) {
            best = total;
            for (int i = 0; i < 3; i++) {
                match[i] = perm[p][i];
                angles[i] = axis_angle_deg(shape[i], inertia[perm[p][i]]);
            }
        }
    }
}

/* ---------- report ---------- */

static void print_matrix(const char *name, const double A[3][3], const char *indent) {
    printf("%s\"%s\": [[%.17g, %.17g, %.17g], [%.17g, %.17g, %.17g], [%.17g, %.17g, %.17g]]",
           indent, name,
           A[0][0], A[0][1], A[0][2],
           A[1][0], A[1][1], A[1][2],
           A[2][0], A[2][1], A[2][2]);
}

static void print_vector(const char *name, const double v[3], const char *indent) {
    printf("%s\"%s\": [%.17g, %.17g, %.17g]", indent, name, v[0], v[1], v[2]);
}

static void report_json(const Mesh *m, const Moments *mom) {
    double lam_shape[3], shape[3][3];
    double lam_I[3], axes[3][3];
    jacobi3(mom->S, lam_shape, shape);
    sort_eigenpairs(lam_shape, shape, 0);      /* long axis first */
    jacobi3(mom->I, lam_I, axes);
    sort_eigenpairs(lam_I, axes, 1);           /* smallest inertia first */

    double lam_origin[3], axes_origin[3][3];
    jacobi3(mom->I_origin, lam_origin, axes_origin);
    sort_eigenpairs(lam_origin, axes_origin, 1);

    double ang[3], ang_origin[3];
    int match[3], match_origin[3];
    match_axes(shape, axes, ang, match);
    match_axes(shape, axes_origin, ang_origin, match_origin);
    double ang_max = fmax(ang[0], fmax(ang[1], ang[2]));
    double ang_origin_max = fmax(ang_origin[0], fmax(ang_origin[1], ang_origin[2]));

    printf("{\n");
    printf("  \"n_vertices\": %d,\n  \"n_triangles\": %d,\n  \"closed\": %s,\n",
           m->nv, m->nf, mom->closed ? "true" : "false");
    printf("  \"volume\": %.17g,\n", mom->V);
    print_vector("centroid", (double[3]){mom->centroid.x, mom->centroid.y, mom->centroid.z}, "  ");
    printf(",\n  \"centroid_offset\": %.17g,\n", norm(mom->centroid));
    print_matrix("second_moment", mom->S, "  ");
    printf(",\n");
    print_matrix("inertia_com", mom->I, "  ");
    printf(",\n");
    print_matrix("inertia_origin", mom->I_origin, "  ");
    printf(",\n");
    printf("  \"parallel_axis_trace\": %.17g,\n", 2.0 * mom->V * dot(mom->centroid, mom->centroid));
    print_vector("principal_moments", lam_I, "  ");
    printf(",\n");
    print_vector("shape_moments", lam_shape, "  ");
    printf(",\n");
    print_matrix("principal_axes", axes, "  ");
    printf(",\n");
    print_matrix("shape_axes", shape, "  ");
    printf(",\n");
    print_vector("axis_angles_deg", ang, "  ");
    printf(",\n  \"axis_angle_max_deg\": %.17g,\n", ang_max);
    print_vector("origin_moments", lam_origin, "  ");
    printf(",\n");
    print_matrix("origin_axes", axes_origin, "  ");
    printf(",\n");
    print_vector("origin_axis_angles_deg", ang_origin, "  ");
    printf(",\n  \"origin_axis_angle_max_deg\": %.17g,\n", ang_origin_max);
    print_vector("axis_match", (double[3]){match[0], match[1], match[2]}, "  ");
    printf("\n}\n");
}

/* ---------- self-test ---------- */

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int near(double a, double b, double tol) { return fabs(a - b) <= tol; }

static void test_cube(void) {
    Mesh m; mesh_init(&m);
    mesh_add_box(&m, -1, -1, -1, 1, 1, 1);      /* 2x2x2, centred */
    Moments mom;
    mesh_moments(&m, &mom);
    check("cube_closed", mom.closed);
    check("cube_volume", near(mom.V, 8.0, 1e-12));
    check("cube_centroid", norm(mom.centroid) < 1e-12);
    /* I = m(dy^2+dz^2)/12 = 8*(4+4)/12 = 16/3 */
    check("cube_inertia_diag", near(mom.I[0][0], 16.0 / 3.0, 1e-12) &&
                               near(mom.I[1][1], 16.0 / 3.0, 1e-12) &&
                               near(mom.I[2][2], 16.0 / 3.0, 1e-12));
    check("cube_offdiag_zero", near(mom.I[0][1], 0.0, 1e-12));
    mesh_free(&m);
}

static void test_offset_box(void) {
    Mesh m; mesh_init(&m);
    mesh_add_box(&m, 9, -2, 3, 13, 2, 9);       /* 4x4x6 box, centroid (11,0,6) */
    Moments mom;
    mesh_moments(&m, &mom);
    check("offset_volume", near(mom.V, 96.0, 1e-12));
    check("offset_centroid", near(mom.centroid.x, 11.0, 1e-12) &&
                             near(mom.centroid.y, 0.0, 1e-12) &&
                             near(mom.centroid.z, 6.0, 1e-12));
    double m_ = 96.0;
    double ix = m_ * (4 * 4 + 6 * 6) / 12.0;
    double iy = m_ * (4 * 4 + 6 * 6) / 12.0;
    double iz = m_ * (4 * 4 + 4 * 4) / 12.0;
    check("offset_inertia_shifted", near(mom.I[0][0], ix, 1e-9) &&
                                    near(mom.I[1][1], iy, 1e-9) &&
                                    near(mom.I[2][2], iz, 1e-9));
    /* the un-shifted tensor must differ: volume is right either way */
    check("offset_origin_tensor_differs",
          !near(mom.I_origin[0][0], mom.I[0][0], 1.0));
    /* trace(I_origin) - trace(I_com) = 2 m |c|^2, c=(11,0,6) -> 314 */
    double tr_origin = mom.I_origin[0][0] + mom.I_origin[1][1] + mom.I_origin[2][2];
    double tr_com = mom.I[0][0] + mom.I[1][1] + mom.I[2][2];
    check("offset_trace_gap", near(tr_origin - tr_com, 2.0 * mom.V * 157.0, 1e-8));
    mesh_free(&m);
}

static void test_rotation(void) {
    Mesh a; mesh_init(&a);
    mesh_add_box(&a, 0, 0, 0, 6, 2, 1);
    Moments ma;
    mesh_moments(&a, &ma);

    /* rotate 90 deg about z: (x,y,z) -> (-y,x,z) */
    Mesh b; mesh_init(&b);
    for (int i = 0; i < a.nv; i++)
        mesh_add_vertex(&b, v3(-a.v[i].y, a.v[i].x, a.v[i].z));
    for (int i = 0; i < a.nf; i++)
        mesh_add_tri(&b, a.f[i].a, a.f[i].b, a.f[i].c);
    Moments mb;
    mesh_moments(&b, &mb);

    double la[3], va[3][3], lb[3], vb[3][3];
    jacobi3(ma.I, la, va); sort_eigenpairs(la, va, 1);
    jacobi3(mb.I, lb, vb); sort_eigenpairs(lb, vb, 1);
    check("rotate_volume", near(ma.V, mb.V, 1e-12));
    /* box centred at (3,1,0.5) maps to (-1,3,0.5) */
    check("rotate_centroid", near(mb.centroid.x, -1.0, 1e-12) &&
                             near(mb.centroid.y, 3.0, 1e-12) &&
                             near(mb.centroid.z, 0.5, 1e-12));
    check("rotate_principal_moments", near(la[0], lb[0], 1e-9) &&
                                      near(la[1], lb[1], 1e-9) &&
                                      near(la[2], lb[2], 1e-9));
    mesh_free(&a); mesh_free(&b);
}

static void test_jacobi(void) {
    /* known eigenpairs: Q has rows as eigenvectors, values 1,2,3 */
    double A[3][3];
    {
        double Q[3][3] = {{1, 0, 0}, {0, 0.6, 0.8}, {0, -0.8, 0.6}};
        double D[3] = {1.0, 2.0, 3.0};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                A[i][j] = 0.0;
                for (int k = 0; k < 3; k++) A[i][j] += Q[i][k] * D[k] * Q[j][k];
            }
    }
    double eval[3], evec[3][3];
    jacobi3(A, eval, evec);
    sort_eigenpairs(eval, evec, 1);
    check("jacobi_values", near(eval[0], 1.0, 1e-12) &&
                           near(eval[1], 2.0, 1e-12) &&
                           near(eval[2], 3.0, 1e-12));
    for (int i = 0; i < 3; i++) {
        double Av[3] = {0, 0, 0};
        for (int r = 0; r < 3; r++)
            for (int s = 0; s < 3; s++) Av[r] += A[r][s] * evec[i][s];
        double res = 0.0;
        for (int r = 0; r < 3; r++) res += (Av[r] - eval[i] * evec[i][r]) *
                                          (Av[r] - eval[i] * evec[i][r]);
        char name[64];
        snprintf(name, sizeof(name), "jacobi_residual_%d", i);
        check(name, sqrt(res) < 1e-12);
    }
}

static void test_open_mesh(void) {
    Mesh m; mesh_init(&m);
    mesh_add_vertex(&m, v3(0, 0, 0));
    mesh_add_vertex(&m, v3(1, 0, 0));
    mesh_add_vertex(&m, v3(0, 1, 0));
    mesh_add_tri(&m, 0, 1, 2);
    Moments mom;
    mesh_moments(&m, &mom);
    check("open_mesh_not_closed", !mom.closed);
    mesh_free(&m);
}

int main(int argc, char **argv) {
    if (argc == 2) {
        Mesh m; mesh_init(&m);
        if (parse_obj(argv[1], &m) != 0) {
            fprintf(stderr, "inertia: cannot read %s\n", argv[1]);
            mesh_free(&m);
            return 2;
        }
        if (m.nf == 0) {
            fprintf(stderr, "inertia: %s has no faces\n", argv[1]);
            mesh_free(&m);
            return 2;
        }
        Moments mom;
        mesh_moments(&m, &mom);
        report_json(&m, &mom);
        mesh_free(&m);
        return mom.closed ? 0 : 3;
    }

    test_cube();
    test_offset_box();
    test_rotation();
    test_jacobi();
    test_open_mesh();

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
