
/* winding: two-sided winding validation for an export -> re-import round trip.
 *
 * The exporter reorders its vertex buffer and canonicalises every index
 * triple to ascending order, which reverses winding for roughly half the
 * faces (parity of the index permutation).  Each boundary runs the same
 * winding_check():
 *
 *   1. shared-edge traversal makes face winding consistent (BFS flips);
 *   2. an area-weighted vote picks the global facing sign.  Each face casts
 *      its area for/against the reference direction r, where r is the +world
 *      axis matching the largest extent of the bounding box.
 *
 * Vote margin:  margin = |W+ - W-| / (W+ + W-)
 * Faces with |n . r| < AXIS_BAND abstain.  margin < TIE_GAP is "close tie":
 * the bbox axis sits in a symmetry-related direction and the normal
 * distribution's dyadic dispersion is reported alongside it.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAXV 320
#define MAXF 640
#define MAXE (3 * MAXF)

#define AXIS_BAND 0.05     /* faces more perpendicular than this abstain */
#define TIE_GAP   0.05     /* margin below this counts as a close tie */
#define DOT_OK    0.99     /* normal agreement threshold (~8 degrees) */

typedef struct { double x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { V3 v[MAXV]; int nv; Tri f[MAXF]; int nf; } Mesh;

/* ---------------- vector helpers ---------------- */

static V3 vsub(V3 p, V3 q) { return (V3){ p.x-q.x, p.y-q.y, p.z-q.z }; }
static V3 vcross(V3 p, V3 q) {
    return (V3){ p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x };
}
static double vdot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static double vlen(V3 p) { return sqrt(vdot(p, p)); }

/* un-normalised face normal; its length is twice the face area */
static V3 face_cross(const Mesh* m, int fi) {
    Tri t = m->f[fi];
    return vcross(vsub(m->v[t.b], m->v[t.a]), vsub(m->v[t.c], m->v[t.a]));
}
static V3 face_normal(const Mesh* m, int fi) {
    V3 c = face_cross(m, fi);
    double l = vlen(c);
    if (l <= 0.0) return (V3){0,0,0};
    return (V3){ c.x/l, c.y/l, c.z/l };
}

/* ---------------- mesh construction ---------------- */

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

static int add_v(Mesh* m, double x, double y, double z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv] = (V3){x, y, z};
    return m->nv++;
}
static int add_f(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a<0||b<0||c<0||a>=m->nv||b>=m->nv||c>=m->nv) return -1;
    if (a==b||b==c||a==c) return -1;
    m->f[m->nf++] = (Tri){a,b,c};
    return m->nf-1;
}

/* Build a surface of revolution around X.  profile rows are (x, radius);
 * a zero radius becomes one shared pole vertex. */
static void build_lathe(Mesh* m, const double prof[][2], int np, int seg) {
    int start[40], len[40];
    for (int i = 0; i < np; i++) {
        double x = prof[i][0], r = prof[i][1];
        if (r <= 1e-12) {
            start[i] = add_v(m, x, 0.0, 0.0);
            len[i] = 1;
        } else {
            start[i] = m->nv;
            len[i] = seg;
            for (int j = 0; j < seg; j++) {
                double th = 2.0 * M_PI * j / seg;
                add_v(m, x, r*cos(th), r*sin(th));
            }
        }
    }
    for (int i = 0; i+1 < np; i++) {
        for (int j = 0; j < seg; j++) {
            int j2 = (j+1) % seg;
            if (len[i] == 1 && len[i+1] > 1) {
                int p = start[i];
                add_f(m, p, start[i+1]+j2, start[i+1]+j);
            } else if (len[i] > 1 && len[i+1] == 1) {
                int p = start[i+1];
                add_f(m, start[i]+j, start[i]+j2, p);
            } else {
                add_f(m, start[i]+j,  start[i]+j2, start[i+1]+j);
                add_f(m, start[i]+j2, start[i+1]+j2, start[i+1]+j);
            }
        }
    }
}

/* Construction-only outward oracle for X-axis lathes: an interior reference
 * point is the axis point at the face centroid's x.  Flips any face whose
 * normal points toward the axis.  Returns number of faces flipped. */
static int orient_outward_radial(Mesh* m) {
    int flipped = 0;
    for (int i = 0; i < m->nf; i++) {
        V3 c = face_cross(m, i);
        Tri t = m->f[i];
        V3 cen = {
            (m->v[t.a].x+m->v[t.b].x+m->v[t.c].x)/3.0,
            (m->v[t.a].y+m->v[t.b].y+m->v[t.c].y)/3.0,
            (m->v[t.a].z+m->v[t.b].z+m->v[t.c].z)/3.0
        };
        V3 axis_pt = { cen.x, 0.0, 0.0 };
        if (vdot(c, vsub(cen, axis_pt)) < 0.0) {
            int tmp = m->f[i].b; m->f[i].b = m->f[i].c; m->f[i].c = tmp;
            flipped++;
        }
    }
    return flipped;
}

static double signed_volume(const Mesh* m) {
    double s = 0.0;
    for (int i = 0; i < m->nf; i++) {
        Tri t = m->f[i];
        s += vdot(m->v[t.a], vcross(m->v[t.b], m->v[t.c]));
    }
    return s / 6.0;
}

/* proper cyclic axis permutation (determinant +1): relabel principal axis */
static void rotate_cyclic(Mesh* m, int mode) {
    for (int i = 0; i < m->nv; i++) {
        V3 p = m->v[i];
        if (mode == 1) m->v[i] = (V3){ p.y, p.z, p.x }; /* old x -> new z */
        if (mode == 2) m->v[i] = (V3){ p.z, p.x, p.y }; /* old x -> new y */
    }
}

/* ---------------- the three models ---------------- */

static void model_rocket(Mesh* m) {
    static const double p[][2] = {
        {-1.60, 0.00}, {-1.52, 0.66}, {-1.25, 0.80},
        { 0.10, 0.80}, { 1.90, 0.42}, { 2.45, 0.00}
    };
    build_lathe(m, p, 6, 20);
    orient_outward_radial(m);
}

static void model_pear(Mesh* m) {
    /* end-for-end mirrored so the canonical +major-axis ray matches outward */
    static const double p[][2] = {
        { 1.00, 0.00}, { 0.88, 0.50}, { 0.50, 0.82}, {-0.05, 0.88},
        {-0.45, 0.72}, {-0.85, 0.40}, {-1.05, 0.12}, {-1.12, 0.00}
    };
    build_lathe(m, p, 8, 20);
    orient_outward_radial(m);
    rotate_cyclic(m, 1);   /* major axis x -> z */
}

/* near-symmetric egg: r(phi) = sin(phi)*(1+k cos(phi)), mild polarity */
static void model_egg(Mesh* m) {
    enum { N = 12 };
    double p[N+1][2];
    const double k = 0.04, sx = 1.12;
    for (int i = 0; i <= N; i++) {
        double phi = M_PI * i / N;
        double r = sin(phi) * (1.0 + k*cos(phi));
        p[i][0] = sx * cos(phi);
        p[i][1] = (i == 0 || i == N) ? 0.0 : r;
    }
    build_lathe(m, p, N+1, 20);
    orient_outward_radial(m);
    rotate_cyclic(m, 2);   /* major axis x -> y */
}

/* ---------------- bbox major axis ---------------- */

typedef struct {
    V3 mn, mx, ext, dir;
    int axis;
} BBox;

static BBox mesh_bbox(const Mesh* m) {
    BBox b;
    b.mn = (V3){ 1e300, 1e300, 1e300 };
    b.mx = (V3){-1e300,-1e300,-1e300 };
    for (int i = 0; i < m->nv; i++) {
        b.mn.x = fmin(b.mn.x, m->v[i].x); b.mx.x = fmax(b.mx.x, m->v[i].x);
        b.mn.y = fmin(b.mn.y, m->v[i].y); b.mx.y = fmax(b.mx.y, m->v[i].y);
        b.mn.z = fmin(b.mn.z, m->v[i].z); b.mx.z = fmax(b.mx.z, m->v[i].z);
    }
    b.ext = vsub(b.mx, b.mn);
    b.axis = 0;
    if (b.ext.y > b.ext.x) b.axis = 1;
    if (b.ext.z > b.ext.x && b.ext.z > b.ext.y) b.axis = 2;
    b.dir = (V3){ b.axis==0, b.axis==1, b.axis==2 };  /* fixed +world sign */
    return b;
}

/* ---------------- symmetric 3x3 eigen-decomposition ---------------- */

static void eigen_sorted(const double Min[3][3], double eval[3], V3 evec[3]) {
    double a[3][3], v[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(a, Min, sizeof(a));
    for (int it = 0; it < 60; it++) {
        double off = 0.0;
        for (int p = 0; p < 3; p++)
            for (int q = p+1; q < 3; q++)
                off += a[p][q]*a[p][q];
        if (off < 1e-22) break;
        for (int p = 0; p < 3; p++) {
            for (int q = p+1; q < 3; q++) {
                double apq = a[p][q];
                if (fabs(apq) < 1e-16) continue;
                double app = a[p][p], aqq = a[q][q];
                double phi = 0.5 * atan2(2.0*apq, aqq-app);
                double c = cos(phi), s = sin(phi);
                for (int k = 0; k < 3; k++) {
                    if (k == p || k == q) continue;
                    double akp = a[k][p], akq = a[k][q];
                    a[k][p] = a[p][k] = c*akp - s*akq;
                    a[k][q] = a[q][k] = s*akp + c*akq;
                }
                a[p][p] = c*c*app - 2*s*c*apq + s*s*aqq;
                a[q][q] = s*s*app + 2*s*c*apq + c*c*aqq;
                a[p][q] = a[q][p] = 0.0;
                for (int k = 0; k < 3; k++) {
                    double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c*vkp - s*vkq;
                    v[k][q] = s*vkp + c*vkq;
                }
            }
        }
    }
    double d[3] = { a[0][0], a[1][1], a[2][2] };
    int idx[3] = {0,1,2};
    for (int i = 0; i < 2; i++)
        for (int j = i+1; j < 3; j++)
            if (d[idx[j]] > d[idx[i]]) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
    for (int i = 0; i < 3; i++) {
        int k = idx[i];
        eval[i] = d[k];
        evec[i] = (V3){ v[0][k], v[1][k], v[2][k] };
    }
}

/* ---------------- area-weighted orientation audit ---------------- */

typedef struct {
    BBox bb;
    double wpos, wneg, wabs, wtot;
    double margin;       /* |W+ - W-| / (W+ + W-) : decisive vote gap */
    double margin_all;   /* |W+ - W-| / total area incl. abstentions */
    int    side;         /* +1 if W+ >= W-, else -1 */
    int    close_tie;
    double cov[3][3];    /* area-weighted dyadic normal dispersion */
    double eval[3];      /* eigenvalues, descending, sum ~ 1 */
    V3     evec[3];      /* principal normal-distribution axes */
    double q_axis;       /* E[(n.r)^2]: 1 axial, 1/3 isotropic, 0 transverse */
    double axis_margin[3]; /* vote margin along each world axis */
} Audit;

static Audit winding_audit(const Mesh* m) {
    Audit z; memset(&z, 0, sizeof(z));
    z.bb = mesh_bbox(m);
    V3 r = z.bb.dir;
    double W = 0.0;
    double ap[3] = {0,0,0}, an[3] = {0,0,0};
    for (int i = 0; i < m->nf; i++) {
        V3 cr = face_cross(m, i);
        double w = vlen(cr);           /* 2*area */
        if (w <= 0.0) continue;
        V3 n = { cr.x/w, cr.y/w, cr.z/w };
        double d = vdot(n, r);
        if (d >  AXIS_BAND) z.wpos += w;
        else if (d < -AXIS_BAND) z.wneg += w;
        else z.wabs += w;
        double nc[3] = { n.x, n.y, n.z };
        for (int ax = 0; ax < 3; ax++) {
            if (nc[ax] >  AXIS_BAND) ap[ax] += w;
            if (nc[ax] < -AXIS_BAND) an[ax] += w;
        }
        W += w;
        for (int a = 0; a < 3; a++) {
            double na = a==0?n.x:a==1?n.y:n.z;
            for (int b2 = 0; b2 < 3; b2++) {
                double nb = b2==0?n.x:b2==1?n.y:n.z;
                z.cov[a][b2] += w * na * nb;
            }
        }
    }
    z.wtot = W;
    if (W > 0.0)
        for (int a = 0; a < 3; a++)
            for (int b2 = 0; b2 < 3; b2++)
                z.cov[a][b2] /= W;
    double decisive = z.wpos + z.wneg;
    z.margin = decisive > 0.0 ? fabs(z.wpos-z.wneg)/decisive : 0.0;
    z.margin_all = W > 0.0 ? fabs(z.wpos-z.wneg)/W : 0.0;
    z.side = z.wpos >= z.wneg ? +1 : -1;
    z.close_tie = z.margin < TIE_GAP;
    z.q_axis = z.cov[z.bb.axis][z.bb.axis];
    for (int ax = 0; ax < 3; ax++) {
        double dec = ap[ax] + an[ax];
        z.axis_margin[ax] = dec > 0.0 ? fabs(ap[ax]-an[ax])/dec : 0.0;
    }
    eigen_sorted(z.cov, z.eval, z.evec);
    return z;
}

/* ---------------- shared-edge winding repair ---------------- */

typedef struct {
    int rel_flips;      /* BFS consistency flips vs as-read order */
    int global_rev;     /* 1 if the axis vote forced a global reversal */
    int check_flips;    /* total faces reversed by the import check */
    int components;
    int conflicts;
    int boundary_edges;
    int final_mask[MAXF];
} RepairInfo;

typedef struct { int lo, hi, fa, ea, fb, eb; } ERec;

/* +1 if face fa's edge ea traverses lo->hi under flip flag f (0 keep,1 rev) */
static int edge_dir_lo_hi(const Mesh* m, const ERec* e, int fa, int ea, int f) {
    const Tri* t = &m->f[fa];
    int idx[3] = { t->a, t->b, t->c };
    int u = idx[ea], v = idx[(ea+1)%3];
    int s = (u == e->lo && v == e->hi) ? +1 : -1;
    return f ? -s : s;
}

static RepairInfo winding_check(Mesh* m) {
    RepairInfo ri; memset(&ri, 0, sizeof(ri));
    static ERec edges[MAXE];
    static int fe[MAXF][3];
    int ne = 0;
    for (int i = 0; i < m->nf; i++)
        for (int k = 0; k < 3; k++)
            fe[i][k] = -1;
    for (int i = 0; i < m->nf; i++) {
        Tri t = m->f[i];
        int idx[3] = { t.a, t.b, t.c };
        for (int k = 0; k < 3; k++) {
            int u = idx[k], v = idx[(k+1)%3];
            int lo = u < v ? u : v, hi = u < v ? v : u;
            int hit = -1;
            for (int e = 0; e < ne; e++)
                if (edges[e].lo == lo && edges[e].hi == hi) { hit = e; break; }
            if (hit < 0) {
                hit = ne++;
                edges[hit] = (ERec){ lo, hi, i, k, -1, -1 };
            } else {
                edges[hit].fb = i; edges[hit].eb = k;
            }
            fe[i][k] = hit;
        }
    }
    for (int e = 0; e < ne; e++)
        if (edges[e].fb < 0) ri.boundary_edges++;

    int assigned[MAXF], flip[MAXF], queue[MAXF];
    memset(assigned, 0, sizeof(assigned));
    memset(flip, 0, sizeof(flip));
    for (int seed = 0; seed < m->nf; seed++) {
        if (assigned[seed]) continue;
        ri.components++;
        int qh = 0, qt = 0;
        assigned[seed] = 1; flip[seed] = 0; queue[qt++] = seed;
        while (qh < qt) {
            int i = queue[qh++];
            for (int k = 0; k < 3; k++) {
                ERec* e = &edges[fe[i][k]];
                int j, kj;
                if (e->fa == i && e->ea == k) { j = e->fb; kj = e->eb; }
                else                            { j = e->fa; kj = e->ea; }
                if (j < 0) continue;
                int di = edge_dir_lo_hi(m, e, i, k, flip[i]);
                /* neighbour must traverse the shared edge oppositely */
                int stored_j = edge_dir_lo_hi(m, e, j, kj, 0);
                int need_fj = (stored_j == -di) ? 0 : 1;
                if (!assigned[j]) {
                    assigned[j] = 1; flip[j] = need_fj; queue[qt++] = j;
                } else if (flip[j] != need_fj) {
                    ri.conflicts++;
                }
            }
        }
    }
    for (int i = 0; i < m->nf; i++) ri.rel_flips += flip[i];

    /* apply relative consistency flips */
    for (int i = 0; i < m->nf; i++) {
        if (!flip[i]) continue;
        int tmp = m->f[i].b; m->f[i].b = m->f[i].c; m->f[i].c = tmp;
    }
    /* area-weighted vote chooses the global sign */
    Audit au = winding_audit(m);
    ri.global_rev = (au.side < 0);
    if (ri.global_rev) {
        for (int i = 0; i < m->nf; i++) {
            int tmp = m->f[i].b; m->f[i].b = m->f[i].c; m->f[i].c = tmp;
        }
    }
    for (int i = 0; i < m->nf; i++) {
        ri.final_mask[i] = flip[i] ^ ri.global_rev;
        ri.check_flips += ri.final_mask[i];
    }
    return ri;
}

/* ---------------- buggy exporter: vertex shuffle + index sort ---------------- */

static void export_reorder(const Mesh* in, Mesh* out, int* truth_mask, unsigned long long seed) {
    mesh_init(out);
    int Q[MAXV], P[MAXV];
    for (int i = 0; i < in->nv; i++) Q[i] = i;
    unsigned long long s = seed;
    for (int i = in->nv-1; i > 0; i--) {     /* Fisher-Yates, fixed LCG */
        s = s*6364136223846793005ULL + 1442695040888963407ULL;
        int j = (int)((s >> 33) % (unsigned long long)(i+1));
        int t = Q[i]; Q[i] = Q[j]; Q[j] = t;
    }
    for (int j = 0; j < in->nv; j++) {
        out->v[j] = in->v[Q[j]];
        P[Q[j]] = j;
    }
    out->nv = in->nv;
    for (int fi = 0; fi < in->nf; fi++) {
        int x[3] = { P[in->f[fi].a], P[in->f[fi].b], P[in->f[fi].c] };
        int inv = 0;
        for (int a = 0; a < 3; a++)
            for (int b = a+1; b < 3; b++)
                if (x[a] > x[b]) inv++;
        truth_mask[fi] = inv & 1;            /* odd permutation reverses winding */
        for (int a = 0; a < 3; a++)
            for (int b = a+1; b < 3; b++)
                if (x[b] < x[a]) { int t = x[a]; x[a] = x[b]; x[b] = t; }
        out->f[fi] = (Tri){ x[0], x[1], x[2] };
    }
    out->nf = in->nf;
}

/* ---------------- minimal OBJ writer / reader ---------------- */

static int write_obj(const char* path, const Mesh* m, const char* name) {
    FILE* fp = fopen(path, "w");
    if (!fp) return -1;
    fprintf(fp, "# winding lab round trip: %s\n", name);
    for (int i = 0; i < m->nv; i++)
        fprintf(fp, "v %.7f %.7f %.7f\n", m->v[i].x, m->v[i].y, m->v[i].z);
    for (int i = 0; i < m->nf; i++)
        fprintf(fp, "f %d %d %d\n", m->f[i].a+1, m->f[i].b+1, m->f[i].c+1);
    fclose(fp);
    return 0;
}

static int read_obj(const char* path, Mesh* m) {
    FILE* fp = fopen(path, "r");
    if (!fp) return -1;
    mesh_init(m);
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == 'v' && (line[1] == ' ' || line[1] == '\t')) {
            double x, y, z;
            if (sscanf(line+1, "%lf %lf %lf", &x, &y, &z) == 3)
                add_v(m, x, y, z);
        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            char* tok = strtok(line+1, " \t\r\n");
            int idx[3], n = 0;
            while (tok && n < 3) {
                idx[n++] = (int)strtol(tok, NULL, 10) - 1; /* strips /vn/vt */
                tok = strtok(NULL, " \t\r\n");
            }
            if (n == 3) add_f(m, idx[0], idx[1], idx[2]);
        }
    }
    fclose(fp);
    return 0;
}

/* ---------------- normal agreement ---------------- */

static void vertex_normals(const Mesh* m, V3* out) {
    for (int i = 0; i < m->nv; i++) out[i] = (V3){0,0,0};
    for (int i = 0; i < m->nf; i++) {
        V3 n = face_cross(m, i);           /* area weighted, unnormalised */
        Tri t = m->f[i];
        out[t.a].x += n.x; out[t.a].y += n.y; out[t.a].z += n.z;
        out[t.b].x += n.x; out[t.b].y += n.y; out[t.b].z += n.z;
        out[t.c].x += n.x; out[t.c].y += n.y; out[t.c].z += n.z;
    }
    for (int i = 0; i < m->nv; i++) {
        double l = vlen(out[i]);
        if (l > 0.0) { out[i].x/=l; out[i].y/=l; out[i].z/=l; }
    }
}

static int nearest_vertex(const Mesh* m, V3 p, double* best_d2) {
    int bi = -1; double bd = 1e300;
    for (int i = 0; i < m->nv; i++) {
        V3 d = vsub(m->v[i], p);
        double d2 = vdot(d, d);
        if (d2 < bd) { bd = d2; bi = i; }
    }
    if (best_d2) *best_d2 = bd;
    return bi;
}

typedef struct { double rate, mean_dot, max_match_err; int agreed; } Agr;

static Agr vertex_agreement(const Mesh* ref, const Mesh* got) {
    V3 nr[MAXV], ng[MAXV];
    vertex_normals(ref, nr);
    vertex_normals(got, ng);
    int ok = 0; double sum = 0.0, maxerr = 0.0;
    for (int i = 0; i < ref->nv; i++) {
        double d2;
        int j = nearest_vertex(got, ref->v[i], &d2);
        if (d2 > maxerr) maxerr = d2;
        double d = (j >= 0) ? vdot(nr[i], ng[j]) : -2.0;
        sum += d;
        if (d >= DOT_OK) ok++;
    }
    Agr a;
    a.rate = ref->nv ? (double)ok/ref->nv : 0.0;
    a.mean_dot = ref->nv ? sum/ref->nv : 0.0;
    a.max_match_err = sqrt(maxerr);
    a.agreed = ok;
    return a;
}

/* face order is preserved by export/import, so compare slot-wise */
static Agr face_agreement_ordered(const Mesh* ref, const Mesh* got) {
    int n = ref->nf < got->nf ? ref->nf : got->nf;
    int ok = 0; double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double d = vdot(face_normal(ref, i), face_normal(got, i));
        sum += d;
        if (d >= DOT_OK) ok++;
    }
    Agr a;
    a.rate = n ? (double)ok/n : 0.0;
    a.mean_dot = n ? sum/n : 0.0;
    a.max_match_err = 0.0;
    a.agreed = ok;
    return a;
}

/* ---------------- test harness ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

typedef struct {
    const char* name;
    int nf, nv, axis;
    double margin;
    int truth_flips, check_flips;
    double raw_face_rate, face_rate, vtx_rate, vtx_mean;
    Audit au;
} Row;

static Row run_case(const char* name, void (*build)(Mesh*),
                    const char* outdir, FILE* tsv) {
    Row row; memset(&row, 0, sizeof(row));
    row.name = name;

    Mesh m0; mesh_init(&m0);
    build(&m0);
    row.nf = m0.nf; row.nv = m0.nv;

    char key[64];
    snprintf(key, sizeof key, "%s_built_limits", name);
    check(key, m0.nv < MAXV && m0.nf < MAXF);

    snprintf(key, sizeof key, "%s_outward_volume", name);
    check(key, signed_volume(&m0) > 0.0);

    /* export-side check on the source: already consistent + vote agrees */
    Mesh me = m0;
    RepairInfo rexp = winding_check(&me);
    snprintf(key, sizeof key, "%s_export_check_zero_flips", name);
    check(key, rexp.check_flips == 0 && rexp.conflicts == 0);

    Audit a0 = winding_audit(&m0);
    row.au = a0; row.margin = a0.margin; row.axis = a0.bb.axis;

    /* buggy export + real file round trip */
    Mesh ms; int truth[MAXF];
    export_reorder(&m0, &ms, truth, 0x9e3779b97f4a7c15ULL ^ (unsigned long long)m0.nv);
    for (int i = 0; i < m0.nf; i++) row.truth_flips += truth[i];

    char path[512];
    snprintf(path, sizeof path, "%s/%s.obj", outdir, name);
    snprintf(key, sizeof key, "%s_obj_written", name);
    check(key, write_obj(path, &ms, name) == 0);

    Mesh mr;
    snprintf(key, sizeof key, "%s_obj_read", name);
    check(key, read_obj(path, &mr) == 0 && mr.nv == ms.nv && mr.nf == ms.nf);

    /* raw re-import: roughly half the face normals are reversed */
    Agr rawf = face_agreement_ordered(&m0, &mr);
    row.raw_face_rate = rawf.rate;

    /* import-side check */
    RepairInfo ri = winding_check(&mr);
    row.check_flips = ri.check_flips;

    int mask_match = (ri.check_flips == row.truth_flips);
    for (int i = 0; mask_match && i < m0.nf; i++)
        if (ri.final_mask[i] != truth[i]) mask_match = 0;
    snprintf(key, sizeof key, "%s_flip_mask_matches_exporter", name);
    check(key, mask_match);
    snprintf(key, sizeof key, "%s_manifold_no_conflict", name);
    check(key, ri.conflicts == 0 && ri.boundary_edges == 0 && ri.components == 1);

    Audit a1 = winding_audit(&mr);
    snprintf(key, sizeof key, "%s_margin_survives_roundtrip", name);
    check(key, fabs(a1.margin - a0.margin) < 1e-6 && a1.side == a0.side);

    Agr ff = face_agreement_ordered(&m0, &mr);
    Agr vv = vertex_agreement(&m0, &mr);
    row.face_rate = ff.rate; row.vtx_rate = vv.rate; row.vtx_mean = vv.mean_dot;

    snprintf(key, sizeof key, "%s_face_normals_restored", name);
    check(key, ff.rate == 1.0);
    snprintf(key, sizeof key, "%s_vertex_match_unique", name);
    check(key, vv.max_match_err < 1e-5);
    snprintf(key, sizeof key, "%s_vertex_normals_restored", name);
    check(key, vv.rate == 1.0);

    double frac = (double)row.truth_flips / m0.nf;
    snprintf(key, sizeof key, "%s_about_half_scrambled", name);
    check(key, frac > 0.35 && frac < 0.65);

    /* detailed record */
    printf("[%s] faces=%d verts=%d axis=%c ext=(%.3f,%.3f,%.3f)\n",
           name, m0.nf, m0.nv, 'x'+a0.bb.axis,
           a0.bb.ext.x, a0.bb.ext.y, a0.bb.ext.z);
    printf("  vote W+=%.1f W-=%.1f abstain=%.1f (area units)  margin=%.4f margin_all=%.4f close_tie=%d\n",
           a0.wpos, a0.wneg, a0.wabs, a0.margin, a0.margin_all, a0.close_tie);
    printf("  exporter_flips=%d/%d (%.3f)  import_check_flips=%d  raw_face_agreement=%.4f\n",
           row.truth_flips, m0.nf, frac, row.check_flips, rawf.rate);
    printf("  corrected face_agreement=%.4f  vertex_agreement=%.4f (mean dot %.6f)\n",
           ff.rate, vv.rate, vv.mean_dot);
    printf("  dispersion eval=(%.4f %.4f %.4f) q_axis=%.4f principal_normal_axis=(%.3f,%.3f,%.3f)\n",
           a0.eval[0], a0.eval[1], a0.eval[2], a0.q_axis,
           a0.evec[0].x, a0.evec[0].y, a0.evec[0].z);
    printf("  axis margins: x=%.6f y=%.6f z=%.6f (bbox major=%c)\n",
           a0.axis_margin[0], a0.axis_margin[1], a0.axis_margin[2], 'x'+a0.bb.axis);
    if (a0.close_tie) {
        int ax = a0.bb.axis;
        double ext_a = ax==0 ? a0.bb.ext.x : ax==1 ? a0.bb.ext.y : a0.bb.ext.z;
        printf("  TIE_DIAG major_axis=%c dir=(%.1f,%.1f,%.1f) extent=%.3f | "
               "vote_margin=%.4f alt_axis_margins=(%.2e,%.2e) | "
               "normal_dispersion eval=(%.4f,%.4f,%.4f) "
               "principal_normal_axis=(%.3f,%.3f,%.3f) q_axis=%.4f\n",
               'x'+ax, a0.bb.dir.x, a0.bb.dir.y, a0.bb.dir.z, ext_a,
               a0.margin,
               a0.axis_margin[(ax+1)%3], a0.axis_margin[(ax+2)%3],
               a0.eval[0], a0.eval[1], a0.eval[2],
               a0.evec[0].x, a0.evec[0].y, a0.evec[0].z, a0.q_axis);
        printf("  TIE-WATCH: reference axis is symmetry-balanced; principal axis and "
               "normal directional dispersion are listed side by side "
               "(isotropic dispersion -> eval ~ 1/3; the global sign is not vote-safe).\n");
    }

    if (tsv) {
        fprintf(tsv, "%s\t%.4f\t%d\t%.4f\n",
                name, row.margin, row.check_flips, row.vtx_rate);
    }
    return row;
}

int main(int argc, char** argv) {
    const char* outdir = "roundtrip_out";
    for (int i = 1; i < argc; i++)
        if (!strcmp(argv[i], "--out") && i+1 < argc) outdir = argv[++i];

    mkdir(outdir, 0755);   /* EEXIST is fine */

    /* tiny math sanity, matching the other labs' self-test style */
    Mesh tri; mesh_init(&tri);
    add_v(&tri,0,0,0); add_v(&tri,1,0,0); add_v(&tri,0,1,0);
    add_f(&tri,0,1,2);
    check("tri_normal_plus_z", fabs(face_normal(&tri,0).z - 1.0) < 1e-12);

    char tsvpath[512];
    snprintf(tsvpath, sizeof tsvpath, "%s/winding_nine.tsv", outdir);
    FILE* tsv = fopen(tsvpath, "w");
    if (tsv) fprintf(tsv, "model\tmargin_ratio\tflipped_faces\tnormal_agreement\n");

    Row r1 = run_case("rocket", model_rocket, outdir, tsv);
    Row r2 = run_case("pear",   model_pear,   outdir, tsv);
    Row r3 = run_case("egg",    model_egg,    outdir, tsv);

    /* the design intent: decisive / medium / close-tie */
    check("rocket_decisive", r1.margin > 0.30 && !r1.au.close_tie);
    check("pear_moderate",   r2.margin > 0.05 && r2.margin < 0.30 && !r2.au.close_tie);
    check("egg_close_tie",   r3.au.close_tie);
    check("egg_dispersion_isotropic",
          r3.au.eval[0] - r3.au.eval[2] < 0.12 &&
          fabs(r3.au.eval[0]+r3.au.eval[1]+r3.au.eval[2]-1.0) < 1e-6);

    if (tsv) fclose(tsv);

    char jsonpath[512];
    snprintf(jsonpath, sizeof jsonpath, "%s/winding_metrics.json", outdir);
    FILE* js = fopen(jsonpath, "w");
    if (js) {
        const Row rows[3] = { r1, r2, r3 };
        fprintf(js, "{\n  \"tie_gap_threshold\": %.2f,\n  \"models\": [\n", TIE_GAP);
        for (int i = 0; i < 3; i++) {
            const Row* r = &rows[i];
            const Audit* a = &r->au;
            int ax = a->bb.axis;
            fprintf(js,
                "    {\"model\":\"%s\",\"faces\":%d,\"vertices\":%d,"
                "\"bbox_major_axis\":\"%c\",\"margin_ratio\":%.6f,"
                "\"margin_ratio_all_area\":%.6f,\"flipped_faces\":%d,"
                "\"exporter_flipped_faces\":%d,\"raw_face_agreement\":%.6f,"
                "\"normal_agreement\":%.6f,\"mean_normal_dot\":%.8f,"
                "\"close_tie\":%s,\"axis_margins\":[%.3e,%.3e,%.3e],"
                "\"normal_dispersion_eval\":[%.6f,%.6f,%.6f],"
                "\"principal_normal_axis\":[%.6f,%.6f,%.6f],\"q_axis\":%.6f}%s\n",
                r->name, r->nf, r->nv, 'x'+ax, r->margin, a->margin_all,
                r->check_flips, r->truth_flips, r->raw_face_rate,
                r->vtx_rate, r->vtx_mean, a->close_tie ? "true" : "false",
                a->axis_margin[0], a->axis_margin[1], a->axis_margin[2],
                a->eval[0], a->eval[1], a->eval[2],
                a->evec[0].x, a->evec[0].y, a->evec[0].z, a->q_axis,
                i < 2 ? "," : "");
        }
        fprintf(js, "  ]\n}\n");
        fclose(js);
    }

    printf("\nnine numbers (3 models x [margin, flips, agreement]):\n");
    printf("model   margin  flips  normal_agreement\n");
    printf("rocket  %.4f  %5d  %.4f\n", r1.margin, r1.check_flips, r1.vtx_rate);
    printf("pear    %.4f  %5d  %.4f\n", r2.margin, r2.check_flips, r2.vtx_rate);
    printf("egg     %.4f  %5d  %.4f%s\n", r3.margin, r3.check_flips, r3.vtx_rate,
           r3.au.close_tie ? "   (close tie; see dispersion)" : "");
    printf("tsv=%s\n", tsvpath);
    printf("json=%s\n", jsonpath);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
