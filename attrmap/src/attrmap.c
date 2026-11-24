
/* attrmap: transfer per-vertex attributes (colors) onto a retopologized
   mesh using closest-point-on-TRIANGLE searches plus barycentric blends.

   The nearest primitive is a face, never just a vertex. Per-face search
   tolerance is set from patch size and target error together:

       h(T)   = longest edge of triangle T           (patch size)
       tol(T) = AM_KH * h(T) + eps_target
       AM_KH  = 1/sqrt(3), vertex covering-radius bound for a triangle

   Two acceptance criteria are measured on three resolution tiers:
   1) round-trip attribute difference at every SOURCE vertex against the
      barycentric sample of the transferred target mesh (max / mean).
   2) degenerate-share: fraction of correspondences whose barycentric max
      weight is >= 1 - AM_ETA, i.e. the "face" mapping collapsed onto one
      vertex and effectively became nearest-vertex copying. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define AM_MAXV 4096
#define AM_MAXF 8192

#define AM_PI 3.14159265358979323846
#define AM_KH (1.0 / 1.73205080756887729353)
#define AM_ETA 0.10          /* barycentric collapse band */
#define AM_EPS_TARGET 1e-3   /* target geometric error, model units */

typedef struct { double x, y, z; } V3;
typedef struct { double r, g, b; } Attr;
typedef struct { int a, b, c; } Tri;

typedef struct {
    V3 v[AM_MAXV]; int nv;
    Tri f[AM_MAXF]; int nf;
    Attr a[AM_MAXV];
} Mesh;

typedef struct {
    V3 p;          /* closest point */
    double w[3];   /* barycentric weights for (a,b,c) of face */
    double dist;   /* absolute distance */
    int face;      /* chosen face index, -1 if mesh has no faces */
} Hit;

static V3 vsub(V3 p, V3 q) { return (V3){ p.x-q.x, p.y-q.y, p.z-q.z }; }
static V3 vadd(V3 p, V3 q) { return (V3){ p.x+q.x, p.y+q.y, p.z+q.z }; }
static V3 vscale(V3 p, double s) { return (V3){ p.x*s, p.y*s, p.z*s }; }
static double vdot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static double vdist(V3 p, V3 q) { V3 d = vsub(p, q); return sqrt(vdot(d, d)); }

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

static int add_vertex(Mesh* m, V3 p, Attr c) {
    if (m->nv >= AM_MAXV) return -1;
    m->v[m->nv] = p;
    m->a[m->nv] = c;
    return m->nv++;
}

static int add_tri(Mesh* m, int a, int b, int c) {
    if (m->nf >= AM_MAXF) return -1;
    if (a<0 || b<0 || c<0 || a>=m->nv || b>=m->nv || c>=m->nv) return -1;
    if (a==b || b==c || a==c) return -1;
    m->f[m->nf] = (Tri){ a, b, c };
    return m->nf++;
}

/* Closest point of p on triangle (a,b,c), Ericson region decomposition.
   Returns the point; w receives the matching barycentric weights. */
static V3 closest_on_tri(V3 p, V3 a, V3 b, V3 c, double w[3]) {
    V3 ab = vsub(b, a), ac = vsub(c, a), ap = vsub(p, a);
    double d1 = vdot(ab, ap), d2 = vdot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        w[0]=1; w[1]=0; w[2]=0; return a;                 /* vertex A */
    }
    V3 bp = vsub(p, b);
    double d3 = vdot(ab, bp), d4 = vdot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        w[0]=0; w[1]=1; w[2]=0; return b;                 /* vertex B */
    }
    V3 cp = vsub(p, c);
    double d5 = vdot(ab, cp), d6 = vdot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        w[0]=0; w[1]=0; w[2]=1; return c;                 /* vertex C */
    }
    double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double v = d1 / (d1 - d3);
        w[0]=1-v; w[1]=v; w[2]=0;                         /* edge AB */
        return vadd(a, vscale(ab, v));
    }
    double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double ww = d2 / (d2 - d6);
        w[0]=1-ww; w[1]=0; w[2]=ww;                       /* edge AC */
        return vadd(a, vscale(ac, ww));
    }
    double va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0) {
        double ww = (d4-d3) / ((d4-d3) + (d5-d6));
        w[0]=0; w[1]=1-ww; w[2]=ww;                       /* edge BC */
        return vadd(b, vscale(vsub(c, b), ww));
    }
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom, ww = vc * denom;
    w[0]=1-v-ww; w[1]=v; w[2]=ww;                         /* interior */
    return vadd(a, vadd(vscale(ab, v), vscale(ac, ww)));
}

/* Brute-force FACE-level closest point. Scans every triangle, so a hit in
   a face interior can never be lost to nearest-vertex preselection. */
static int mesh_closest(const Mesh* m, V3 p, Hit* hit) {
    hit->face = -1;
    hit->dist = 1e300;
    if (m->nf == 0) return -1;
    for (int i = 0; i < m->nf; i++) {
        Tri t = m->f[i];
        double w[3];
        V3 q = closest_on_tri(p, m->v[t.a], m->v[t.b], m->v[t.c], w);
        double d = vdist(p, q);
        if (d < hit->dist) {
            hit->dist = d; hit->p = q; hit->face = i;
            hit->w[0]=w[0]; hit->w[1]=w[1]; hit->w[2]=w[2];
        }
    }
    return 0;
}

/* characteristic patch size of one face: longest edge */
static double face_size(const Mesh* m, int fi) {
    Tri t = m->f[fi];
    double e0 = vdist(m->v[t.a], m->v[t.b]);
    double e1 = vdist(m->v[t.b], m->v[t.c]);
    double e2 = vdist(m->v[t.c], m->v[t.a]);
    double h = e0 > e1 ? e0 : e1;
    return h > e2 ? h : e2;
}

/* tol = h/sqrt(3) + eps_target. Returns -1 on bad arguments. */
static double search_tol(double h, double eps_target) {
    if (h < 0.0 || eps_target < 0.0 || !isfinite(h) || !isfinite(eps_target)) return -1.0;
    return AM_KH * h + eps_target;
}

/* A mapping is degenerate when one barycentric weight owns the blend. */
static int hit_is_degenerate(const Hit* h) {
    double wm = h->w[0] > h->w[1] ? h->w[0] : h->w[1];
    if (h->w[2] > wm) wm = h->w[2];
    return wm >= 1.0 - AM_ETA;
}

static Attr attr_blend(Attr ca, Attr cb, Attr cc, const double w[3]) {
    return (Attr){
        w[0]*ca.r + w[1]*cb.r + w[2]*cc.r,
        w[0]*ca.g + w[1]*cb.g + w[2]*cc.g,
        w[0]*ca.b + w[1]*cb.b + w[2]*cc.b
    };
}

static double attr_dist(Attr x, Attr y) {
    double dr = x.r-y.r, dg = x.g-y.g, db = x.b-y.b;
    return sqrt(dr*dr + dg*dg + db*db);
}

/* analytic smooth color field, shared by both LODs */
static Attr shade(V3 p) {
    return (Attr){
        0.5 + 0.5*sin(2.0*AM_PI*p.x),
        0.5 + 0.5*cos(2.0*AM_PI*p.y),
        sin(AM_PI*p.x)*sin(AM_PI*p.y)
    };
}

/* Regular triangulated height-field grid, K cells per side.
   Two triangles per cell with a fixed diagonal. */
static int build_grid(Mesh* m, int K, double amp) {
    mesh_init(m);
    if (K < 1) return -1;
    int stride = K + 1;
    for (int i = 0; i <= K; i++) {
        for (int j = 0; j <= K; j++) {
            double x = (double)i / K, y = (double)j / K;
            V3 p = { x, y, amp*sin(AM_PI*x)*sin(AM_PI*y) };
            if (add_vertex(m, p, shade(p)) < 0) return -1;
        }
    }
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            int a = i*stride + j;
            int b = (i+1)*stride + j;
            int c = (i+1)*stride + (j+1);
            int d = i*stride + (j+1);
            if (add_tri(m, a, b, d) < 0) return -1;
            if (add_tri(m, b, c, d) < 0) return -1;
        }
    }
    return 0;
}

static double mean_face_size(const Mesh* m) {
    double s = 0.0;
    for (int i = 0; i < m->nf; i++) s += face_size(m, i);
    return m->nf ? s / m->nf : 0.0;
}

typedef struct {
    int S, D;
    int nsrc, nfsrc, ndst, nfdst;
    double face_ratio, vratio;
    double hsrc, hdst, toldst;
    /* criterion 1 + 2, measured over SOURCE vertices on the DST mesh */
    double max_diff, mean_diff, degenerate_pct;
    int rejected;
    /* nearest-vertex-only control on the same population */
    double vmax_diff, vmean_diff, vdegenerate_pct;
} TierResult;

static void run_tier(int S, int D, TierResult* r) {
    Mesh src, dst;
    memset(r, 0, sizeof(*r));
    r->S = S; r->D = D;

    build_grid(&src, S, 0.08);
    build_grid(&dst, D, 0.08);
    r->nsrc = src.nv; r->nfsrc = src.nf;
    r->ndst = dst.nv; r->nfdst = dst.nf;
    r->face_ratio = (double)src.nf / dst.nf;
    r->vratio = (double)src.nv / dst.nv;
    r->hsrc = mean_face_size(&src);
    r->hdst = mean_face_size(&dst);
    r->toldst = search_tol(r->hdst, AM_EPS_TARGET);

    /* Forward transfer: each retopo vertex samples the dense SOURCE faces.
       (Nested grids make these exact shared vertices in this experiment.) */
    for (int i = 0; i < dst.nv; i++) {
        Hit h;
        if (mesh_closest(&src, dst.v[i], &h) != 0) continue;
        double tol = search_tol(face_size(&src, h.face), AM_EPS_TARGET);
        if (h.dist > tol) continue;
        Tri t = src.f[h.face];
        dst.a[i] = attr_blend(src.a[t.a], src.a[t.b], src.a[t.c], h.w);
    }

    /* Criterion 1 + 2: every SOURCE vertex back onto DST faces. */
    double sum = 0.0, mx = 0.0;
    long ndeg = 0;
    for (int i = 0; i < src.nv; i++) {
        Hit h;
        mesh_closest(&dst, src.v[i], &h);
        double tol = search_tol(face_size(&dst, h.face), AM_EPS_TARGET);
        if (h.dist > tol) r->rejected++;
        Tri t = dst.f[h.face];
        Attr got = attr_blend(dst.a[t.a], dst.a[t.b], dst.a[t.c], h.w);
        double e = attr_dist(src.a[i], got);
        sum += e;
        if (e > mx) mx = e;
        if (hit_is_degenerate(&h)) ndeg++;
    }
    r->max_diff = mx;
    r->mean_diff = sum / src.nv;
    r->degenerate_pct = 100.0 * (double)ndeg / src.nv;

    /* Control: vertex-only nearest search, same source population. */
    double vsum = 0.0, vmx = 0.0;
    long vdeg = 0;
    for (int i = 0; i < src.nv; i++) {
        int best = -1; double bd = 1e300;
        for (int j = 0; j < dst.nv; j++) {
            double dd = vdist(src.v[i], dst.v[j]);
            if (dd < bd) { bd = dd; best = j; }
        }
        double e = attr_dist(src.a[i], dst.a[best]);
        vsum += e;
        if (e > vmx) vmx = e;
        vdeg++;   /* pure nearest-vertex copy is degenerate by construction */
    }
    r->vmax_diff = vmx;
    r->vmean_diff = vsum / src.nv;
    r->vdegenerate_pct = 100.0 * (double)vdeg / src.nv;
}

/* ---------------- self-test harness (repo convention) ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int near_eq(double x, double y, double eps) { return fabs(x-y) <= eps; }

int main(void) {
    /* --- unit: closest point regions and barycentric weights --- */
    V3 A={0,0,0}, B={1,0,0}, C={0,1,0};
    double w[3];

    V3 qi = closest_on_tri((V3){0.2,0.2,1.0}, A,B,C,w);
    check("interior_point_xy", near_eq(qi.x,0.2,1e-9) && near_eq(qi.y,0.2,1e-9) && near_eq(qi.z,0.0,1e-9));
    check("interior_bary", near_eq(w[0],0.6,1e-9) && near_eq(w[1],0.2,1e-9) && near_eq(w[2],0.2,1e-9));
    check("interior_not_degenerate", (w[0]>0.5 && w[1]>0.1 && w[2]>0.1));

    V3 qv = closest_on_tri((V3){-1,-1,0}, A,B,C,w);
    check("vertex_region", near_eq(qv.x,0,1e-9) && near_eq(qv.y,0,1e-9) && w[0] > 0.999);

    V3 qe = closest_on_tri((V3){0.6,0.6,0}, A,B,C,w);
    check("edge_region", near_eq(qe.x,0.5,1e-9) && near_eq(qe.y,0.5,1e-9));
    check("edge_bary", near_eq(w[0],0.0,1e-9) && near_eq(w[1],0.5,1e-9) && near_eq(w[2],0.5,1e-9));

    check("bary_sums_one", near_eq(w[0]+w[1]+w[2],1.0,1e-9));

    /* face blend must beat nearest-vertex copy on an interior query */
    Attr red={1,0,0}, green={0,1,0}, blue={0,0,1};
    double wi[3];
    closest_on_tri((V3){0.34,0.33,0}, A,B,C,wi);
    Attr blended = attr_blend(red,green,blue,wi);
    double d_to_red = attr_dist(blended, red);
    check("face_blend_not_vertex_color", d_to_red > 0.5);
    check("face_blend_balanced", wi[0] > 0.2 && wi[1] > 0.2 && wi[2] > 0.2);

    /* --- tolerance from patch size + target error --- */
    double t = search_tol(0.30, 1e-3);
    check("tol_formula", near_eq(t, AM_KH*0.30 + 1e-3, 1e-12));
    check("tol_grows_with_patch", search_tol(0.40,1e-3) > search_tol(0.20,1e-3));
    check("tol_target_floor", search_tol(0.0,2e-3) >= 2e-3 - 1e-12);
    check("tol_bad_args", search_tol(-1.0,1e-3) < 0.0 && search_tol(0.3,-1.0) < 0.0);

    /* --- degeneracy classifier --- */
    Hit hd; memset(&hd, 0, sizeof(hd));
    hd.w[0]=1.0; hd.w[1]=0; hd.w[2]=0;
    check("degen_vertex_copy", hit_is_degenerate(&hd));
    hd.w[0]=0.89; hd.w[1]=0.11; hd.w[2]=0.0;
    check("degen_band_boundary_out", !hit_is_degenerate(&hd));
    hd.w[0]=0.91; hd.w[1]=0.05; hd.w[2]=0.04;
    check("degen_band_boundary_in", hit_is_degenerate(&hd));

    /* --- grid builder counts --- */
    Mesh g;
    build_grid(&g, 8, 0.08);
    check("grid_vertex_count", g.nv == 81);
    check("grid_face_count", g.nf == 128);
    check("mean_face_size_sane", mean_face_size(&g) > 0.15 && mean_face_size(&g) < 0.20);

    /* --- three tiers: face ratios 4x / 16x / 64x --- */
    TierResult r[3];
    run_tier(16, 8, &r[0]);
    run_tier(24, 6, &r[1]);
    run_tier(32, 4, &r[2]);

    for (int k = 0; k < 3; k++) {
        char nm[64];
        snprintf(nm, sizeof(nm), "tier%d_ratio_exact", k+1);
        double expect_ratio[3] = {4.0, 16.0, 64.0};
        check(nm, near_eq(r[k].face_ratio, expect_ratio[k], 1e-9));
        snprintf(nm, sizeof(nm), "tier%d_no_rejected", k+1);
        check(nm, r[k].rejected == 0);
        snprintf(nm, sizeof(nm), "tier%d_errors_finite_positive", k+1);
        check(nm, isfinite(r[k].max_diff) && r[k].max_diff > 0.0 &&
                  r[k].mean_diff > 0.0 && r[k].mean_diff <= r[k].max_diff);
    }

    /* criterion 2 diagnostic: largest-ratio tier must NOT be the most
       degenerate under true face-level search */
    check("crit_maxratio_not_most_degenerate",
          r[2].degenerate_pct < r[0].degenerate_pct &&
          r[2].degenerate_pct < r[1].degenerate_pct);
    check("crit_extreme_ratio_low_degenerate", r[2].degenerate_pct < 10.0);
    check("crit_degen_decreasing",
          r[0].degenerate_pct > r[1].degenerate_pct &&
          r[1].degenerate_pct > r[2].degenerate_pct);

    /* face-level transfer must beat vertex-only copying at extreme ratio */
    check("crit_face_beats_vertex_extreme",
          r[2].mean_diff < r[2].vmean_diff &&
          r[2].max_diff < r[2].vmax_diff);
    check("control_vertex_only_full_degenerate",
          r[0].vdegenerate_pct == 100.0 &&
          r[1].vdegenerate_pct == 100.0 &&
          r[2].vdegenerate_pct == 100.0);

    /* tolerance tracks the growing retopo patch size across tiers */
    check("tol_scales_with_patch_across_tiers",
          r[0].hdst < r[1].hdst && r[1].hdst < r[2].hdst &&
          r[0].toldst < r[1].toldst && r[1].toldst < r[2].toldst);

    /* ---------------- machine-greppable nine numbers ---------------- */
    for (int k = 0; k < 3; k++) {
        printf("tier%d_max_attr_diff=%.6f\n",  k+1, r[k].max_diff);
        printf("tier%d_mean_attr_diff=%.6f\n", k+1, r[k].mean_diff);
        printf("tier%d_degenerate_pct=%.2f\n", k+1, r[k].degenerate_pct);
    }

    /* ---------------- human-readable experiment tables -------------- */
    printf("\n== nine numbers: attribute transfer, face-level mapping ==\n");
    printf("tier  S/D   face_ratio  max_diff  mean_diff  degenerate%%\n");
    for (int k = 0; k < 3; k++) {
        printf("%d     %2d/%-2d %6.0fx    %.4f    %.4f     %6.2f%%\n",
               k+1, r[k].S, r[k].D, r[k].face_ratio,
               r[k].max_diff, r[k].mean_diff, r[k].degenerate_pct);
    }

    printf("\n== patch size vs search tolerance (target eps=%.0e, kh=1/sqrt3) ==\n",
           AM_EPS_TARGET);
    printf("tier  face_ratio  h_src    h_dst    tol_dst\n");
    for (int k = 0; k < 3; k++) {
        printf("%d     %6.0fx    %.4f   %.4f   %.4f\n",
               k+1, r[k].face_ratio, r[k].hsrc, r[k].hdst, r[k].toldst);
    }

    printf("\n== control: vertex-only nearest search (same vertices) ==\n");
    printf("tier  vmax_diff  vmean_diff  degenerate%%\n");
    for (int k = 0; k < 3; k++) {
        printf("%d     %.4f     %.4f      %6.2f%%\n",
               k+1, r[k].vmax_diff, r[k].vmean_diff, r[k].vdegenerate_pct);
    }

    printf("\ndiagnostic max_ratio_highest_degenerate=%d (1 would flag vertex-only search)\n",
           (r[2].degenerate_pct > r[0].degenerate_pct &&
            r[2].degenerate_pct > r[1].degenerate_pct));

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
