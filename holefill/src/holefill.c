
/* holefill: orient a hole fan from the boundary ring's own cross-product
 * sign, then gate the result by added-area vs projected-hole-area.
 *
 * Orientation rule (no world-up assumption):
 *   Boundary edges are taken in the direction they have inside their
 *   incident (outward-wound) triangles and chained into a directed ring L.
 *   For an outward-wound mesh the signed ring normal
 *       n_loop = norm( sum_i cross(p_i, p_{i+1}) )        (Newell)
 *              || sign of sum_i cross(e_i, e_{i+1})       (adjacent edges)
 *   points INWARD, i.e. opposite the patch outward normal.  Each fan
 *   triangle (center, L[i], L[i+1]) is therefore flipped whenever
 *       dot(normal_tri, -n_loop) < 0.
 *   Rotating the whole model rotates every vector here, so the decision is
 *   invariant to hole orientation in the world.
 *
 * Second gate:
 *   r = sum(fill triangle areas) / projected hole area (= |Newell|/2).
 *   A planar centroid fan gives r ~= 1; an apex outside / balloon patch is
 *   rejected when r > AREA_RATIO_MAX.  Orientation and area must both pass.
 *
 * The program also contains the old buggy policy ORIENT_FIXED_UP
 * ("flip if dot(normal, (0,0,1)) < 0") to reproduce the all-black downward
 * hole and the reference-direction diagnostic.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#define MAXV   48
#define MAXF   64
#define LOOPV  16
#define MAXBE  48
#define MAXLOOPS 4

#define RING_R 1.0f
#define OUTER_R 2.0f
#define AREA_RATIO_MAX 1.5f

#define RW 128
#define RH 128
#define PX_PER_UNIT 30.0f
#define ROI_FACTOR 0.80f      /* disk strictly inside the r*cos30 inradius */
#define BLACK_LUM 40.0f

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct {
    V3 v[MAXV]; int nv;
    Tri f[MAXF]; int nf;
} Mesh;

static V3 vsub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 vadd(V3 p, V3 q) { V3 r = { p.x+q.x, p.y+q.y, p.z+q.z }; return r; }
static V3 vscale(V3 p, float s) { V3 r = { p.x*s, p.y*s, p.z*s }; return r; }
static V3 vcross(V3 p, V3 q) {
    V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x };
    return r;
}
static float vdot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float vlen(V3 p) { return sqrtf(vdot(p,p)); }
static V3 vnorm(V3 p) {
    float l = vlen(p);
    V3 z = {0,0,0};
    if (l <= 1e-20f) return z;
    return vscale(p, 1.0f/l);
}

static int add_v(Mesh* m, V3 p) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv] = p;
    return m->nv++;
}
static int add_f(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a<0||b<0||c<0||a>=m->nv||b>=m->nv||c>=m->nv) return -1;
    if (a==b||b==c||a==c) return -1;
    m->f[m->nf].a=a; m->f[m->nf].b=b; m->f[m->nf].c=c;
    return m->nf++;
}
static V3 face_raw_normal(const Mesh* m, int fi) {
    Tri t = m->f[fi];
    V3 a = m->v[t.a], b = m->v[t.b], c = m->v[t.c];
    return vcross(vsub(b,a), vsub(c,a));
}

/* ---------------- oriented test models ---------------- */

typedef enum { HOLE_UP = 0, HOLE_SIDE = 1, HOLE_DOWN = 2 } HoleKind;

typedef struct {
    const char* name;
    V3 u, v, n;       /* orthonormal, u x v = n = true outward normal */
} Basis;

static Basis basis_of(HoleKind k) {
    Basis b;
    if (k == HOLE_UP) {
        b.name = "up";
        b.u = (V3){1,0,0}; b.v = (V3){0,1,0}; b.n = (V3){0,0,1};
    } else if (k == HOLE_SIDE) {
        b.name = "side";
        b.u = (V3){0,1,0}; b.v = (V3){0,0,1}; b.n = (V3){1,0,0};
    } else {
        b.name = "down";
        b.u = (V3){1,0,0}; b.v = (V3){0,-1,0}; b.n = (V3){0,0,-1};
    }
    return b;
}

static V3 ring_point(const Basis* b, float radius, float theta) {
    V3 d = vadd(vscale(b->u, cosf(theta)*radius),
                vscale(b->v, sinf(theta)*radius));
    return d; /* rings are centered on origin */
}

/* Hexagonal hole (inner radius 1) in a hexagonal flange (outer radius 2).
 * All 12 flange tris are wound so their geometric normal is +b.n. */
static Mesh build_flange(HoleKind kind) {
    Mesh m; memset(&m, 0, sizeof(m));
    Basis b = basis_of(kind);
    enum { K = 6 };
    int inner[K], outer[K];
    for (int i = 0; i < K; i++) {
        float th = (float)i * (2.0f*3.14159265358979323846f / K);
        inner[i] = add_v(&m, ring_point(&b, RING_R, th));
        outer[i] = add_v(&m, ring_point(&b, OUTER_R, th));
    }
    for (int i = 0; i < K; i++) {
        int j = (i+1) % K;
        add_f(&m, inner[i], outer[i], outer[j]);   /* normal +n */
        add_f(&m, inner[i], outer[j], inner[j]);   /* inner edge inner[j]->inner[i] */
    }
    return m;
}

/* ---------------- boundary loops ---------------- */

typedef struct { int a, b, owner; } DEdge;
typedef struct { int idx[LOOPV]; int n; } Loop;

/* Edges with no reverse-directed partner: boundary edges, directed exactly
 * as they occur in the owner triangle's winding. */
static int boundary_directed(const Mesh* m, DEdge* out) {
    int ne = 0;
    for (int fi = 0; fi < m->nf; fi++) {
        Tri t = m->f[fi];
        int ea[3] = { t.a, t.b, t.c };
        int eb[3] = { t.b, t.c, t.a };
        for (int e = 0; e < 3; e++) {
            int x = ea[e], y = eb[e];
            int reverse = 0;
            for (int fj = 0; fj < m->nf && !reverse; fj++) {
                Tri q = m->f[fj];
                int qa[3] = { q.a, q.b, q.c };
                int qb[3] = { q.b, q.c, q.a };
                for (int g = 0; g < 3; g++) {
                    if (qa[g] == y && qb[g] == x) { reverse = 1; break; }
                }
            }
            if (!reverse) {
                if (ne < MAXBE) out[ne++] = (DEdge){x, y, fi};
            }
        }
    }
    return ne;
}

static int extract_loops(DEdge* be, int ne, Loop* loops) {
    int used[MAXBE]; memset(used, 0, sizeof(used));
    int nl = 0;
    for (int s = 0; s < ne; s++) {
        if (used[s]) continue;
        if (nl >= MAXLOOPS) break;
        Loop* L = &loops[nl]; memset(L, 0, sizeof(*L));
        int start = be[s].a, cur = start;
        for (;;) {
            if (L->n >= LOOPV) break;
            L->idx[L->n++] = cur;
            int pick = -1;
            for (int j = 0; j < ne; j++) {
                if (!used[j] && be[j].a == cur) { pick = j; break; }
            }
            if (pick < 0) break;
            used[pick] = 1;
            cur = be[pick].b;
            if (cur == start) break;
        }
        nl++;
    }
    return nl;
}

static V3 loop_newell(const Mesh* m, const Loop* L) {
    V3 N = {0,0,0};
    for (int i = 0; i < L->n; i++) {
        V3 p = m->v[L->idx[i]];
        V3 q = m->v[L->idx[(i+1) % L->n]];
        N.x += p.y*q.z - p.z*q.y;
        N.y += p.z*q.x - p.x*q.z;
        N.z += p.x*q.y - p.y*q.x;
    }
    return N; /* magnitude = 2 * projected area */
}

/* sum of cross(e_i, e_{i+1}): same signed direction as Newell for a
 * convex planar ring; this is the literal adjacent-edge cross sign. */
static V3 loop_adjacent_cross(const Mesh* m, const Loop* L) {
    V3 C = {0,0,0};
    for (int i = 0; i < L->n; i++) {
        V3 p0 = m->v[L->idx[i]];
        V3 p1 = m->v[L->idx[(i+1) % L->n]];
        V3 p2 = m->v[L->idx[(i+2) % L->n]];
        V3 e1 = vsub(p1, p0), e2 = vsub(p2, p1);
        V3 x = vcross(e1, e2);
        C = vadd(C, x);
    }
    return C;
}

static float loop_projected_area(const Mesh* m, const Loop* L) {
    return 0.5f * vlen(loop_newell(m, L));
}

static Loop select_hole_loop(const Mesh* m, const Loop* loops, int nl) {
    int best = 0;
    float bestA = 1e30f;
    for (int i = 0; i < nl; i++) {
        float A = loop_projected_area(m, &loops[i]);
        if (A < bestA) { bestA = A; best = i; }   /* inner ring = smaller */
    }
    return loops[best];
}

/* ---------------- hole filling ---------------- */

typedef enum { ORIENT_RING_SIGN = 0, ORIENT_FIXED_UP = 1 } OrientMode;

static const V3 FIXED_UP = {0,0,1};

/* Outward winding target read from the directed ring itself: minus the sign
 * of the adjacent-edge cross sum (== minus Newell for a planar ring). */
static V3 ring_outward_target(const Mesh* m, const Loop* L) {
    V3 c = loop_adjacent_cross(m, L);
    if (vlen(c) < 1e-12f) c = loop_newell(m, L);
    return vscale(vnorm(c), -1.0f);
}

/* Append centroid fan. Returns center vertex index, sets fstart/fcount.
 * ORIENT_RING_SIGN: target normal from the ring cross sign (mesh intrinsic).
 * ORIENT_FIXED_UP : old buggy policy, target = world +Z. */
static int fill_centroid_fan(Mesh* m, const Loop* L, OrientMode mode,
                             int* fstart, int* fcount) {
    V3 c = {0,0,0};
    for (int i = 0; i < L->n; i++) c = vadd(c, m->v[L->idx[i]]);
    c = vscale(c, 1.0f / (float)L->n);
    int ci = add_v(m, c);

    V3 target = (mode == ORIENT_RING_SIGN) ? ring_outward_target(m, L)
                                           : FIXED_UP;

    *fstart = m->nf;
    for (int i = 0; i < L->n; i++) {
        int a = L->idx[i], b = L->idx[(i+1) % L->n];
        /* normal of candidate tri (center, a, b) */
        V3 nt = vcross(vsub(m->v[a], c), vsub(m->v[b], c));
        /* dot<0 -> swap b,a so the tri faces target. dot==0 (ring plane
         * parallel to fixed reference) keeps the candidate ordering. */
        if (vdot(nt, target) < 0.0f) add_f(m, ci, b, a);
        else                          add_f(m, ci, a, b);
    }
    *fcount = m->nf - *fstart;
    return ci;
}

/* Deliberately over-sized fill: apex lifted h along outward normal.
 * Ring-sign still orients every tri outward; area gate must reject it. */
static int fill_balloon(Mesh* m, const Loop* L, V3 n, float h,
                        int* fstart, int* fcount) {
    V3 c = {0,0,0};
    for (int i = 0; i < L->n; i++) c = vadd(c, m->v[L->idx[i]]);
    c = vscale(c, 1.0f/(float)L->n);
    c = vadd(c, vscale(n, h));
    int ci = add_v(m, c);
    V3 target = ring_outward_target(m, L);
    *fstart = m->nf;
    for (int i = 0; i < L->n; i++) {
        int a = L->idx[i], b = L->idx[(i+1) % L->n];
        V3 ta = m->v[a], tb = m->v[b];
        V3 nt = vcross(vsub(ta, c), vsub(tb, c));
        if (vdot(nt, target) < 0.0f) add_f(m, ci, b, a);
        else                          add_f(m, ci, a, b);
    }
    *fcount = m->nf - *fstart;
    return ci;
}

/* ---------------- metrics ---------------- */

typedef struct {
    float dir_ratio;       /* fill tris with normal dot truth > 0, fraction */
    int   dir_ok;
    int   outward_v;        /* patch vertices whose fill-only normal is out */
    int   patch_v;
    int   outward_all;      /* whole-mesh area-weighted normals pointing out */
    int   nv_all;
    int   patch_out_allmesh;/* patch vertices under whole-mesh normals */
    float area_fill;
    float area_proj;
    float area_ratio;
    int   area_ok;
    int   black_px;         /* dark pixels inside the hole ROI */
    int   roi_px;
} Metrics;

static float fill_area(const Mesh* m, int fstart, int fcount) {
    float A = 0.0f;
    for (int f = fstart; f < fstart+fcount; f++)
        A += 0.5f * vlen(face_raw_normal(m, f));
    return A;
}

static float direction_ratio(const Mesh* m, int fstart, int fcount, V3 truth) {
    int good = 0;
    for (int f = fstart; f < fstart+fcount; f++)
        if (vdot(face_raw_normal(m, f), truth) > 0.0f) good++;
    return fcount ? (float)good/(float)fcount : 0.0f;
}

/* area-weighted normals accumulated from fill faces only, counted over the
 * center + ring vertices that make up the patch */
static int outward_patch_vertices(const Mesh* m, const Loop* L, int ci,
                                  int fstart, int fcount, V3 truth) {
    V3 acc[MAXV]; memset(acc, 0, sizeof(acc));
    for (int f = fstart; f < fstart+fcount; f++) {
        Tri t = m->f[f];
        V3 n = face_raw_normal(m, f);
        acc[t.a] = vadd(acc[t.a], n);
        acc[t.b] = vadd(acc[t.b], n);
        acc[t.c] = vadd(acc[t.c], n);
    }
    int out = 0;
    if (vdot(acc[ci], truth) > 0.0f) out++;
    for (int i = 0; i < L->n; i++)
        if (vdot(acc[L->idx[i]], truth) > 0.0f) out++;
    return out;
}

/* whole-mesh area-weighted vertex normals (same convention as normal.c);
 * returns total outward vertices and, separately, outward patch vertices */
static void outward_full_mesh(const Mesh* m, const Loop* L, int ci, V3 truth,
                              int* out_all, int* out_patch) {
    V3 acc[MAXV]; memset(acc, 0, sizeof(acc));
    for (int f = 0; f < m->nf; f++) {
        Tri t = m->f[f];
        V3 n = face_raw_normal(m, f);
        acc[t.a] = vadd(acc[t.a], n);
        acc[t.b] = vadd(acc[t.b], n);
        acc[t.c] = vadd(acc[t.c], n);
    }
    int all = 0;
    for (int i = 0; i < m->nv; i++)
        if (vdot(acc[i], truth) > 0.0f) all++;
    int patch = (vdot(acc[ci], truth) > 0.0f);
    for (int i = 0; i < L->n; i++)
        if (vdot(acc[L->idx[i]], truth) > 0.0f) patch++;
    *out_all = all;
    *out_patch = patch;
}

/* ---------------- tiny orthographic software renderer ---------------- */

static void project(V3 p, V3 g, V3 bu, V3 bv, float* sx, float* sy) {
    V3 d = vsub(p, g);
    *sx = RW*0.5f + PX_PER_UNIT*vdot(d, bu);
    *sy = RH*0.5f - PX_PER_UNIT*vdot(d, bv);   /* screen y grows downward */
}

static float edge2(float ax, float ay, float bx, float by,
                   float px, float py) {
    return (bx-ax)*(py-ay) - (by-ay)*(px-ax);
}

/* head-on orthographic camera on the outward side, headlight tilted a touch.
 * Shading is pure Lambert against a black background: inward normals get
 * dot <= 0 and render exactly black (no ambient, no culling needed). */
static void render_mesh(const Mesh* m, int fstart, int fcount,
                        V3 bu, V3 bv, V3 bn, unsigned char* rgb) {
    memset(rgb, 0, (size_t)RW*RH*3);
    V3 Ld = vnorm(vadd(vadd(bn, vscale(bu, 0.30f)), vscale(bv, 0.20f)));

    for (int fi = 0; fi < m->nf; fi++) {
        Tri t = m->f[fi];
        int is_fill = (fi >= fstart && fi < fstart+fcount);
        const float* base = is_fill ? (const float[]){255,176,64}
                                    : (const float[]){178,184,194};
        V3 fn = vnorm(face_raw_normal(m, fi));
        float diff = vdot(fn, Ld);
        float shade = diff > 0.0f ? diff : 0.0f;   /* inward -> 0 -> black */
        float col[3] = { base[0]*shade, base[1]*shade, base[2]*shade };

        V3 g0 = {0,0,0};
        float x0,y0,x1,y1,x2,y2;
        project(m->v[t.a], g0, bu, bv, &x0,&y0);
        project(m->v[t.b], g0, bu, bv, &x1,&y1);
        project(m->v[t.c], g0, bu, bv, &x2,&y2);

        int ix0 = (int)floorf(fminf(x0,fminf(x1,x2)));
        int ix1 = (int)ceilf (fmaxf(x0,fmaxf(x1,x2)));
        int iy0 = (int)floorf(fminf(y0,fminf(y1,y2)));
        int iy1 = (int)ceilf (fmaxf(y0,fmaxf(y1,y2)));
        if (ix0 < 0) ix0 = 0;
        if (iy0 < 0) iy0 = 0;
        if (ix1 >= RW) ix1 = RW-1;
        if (iy1 >= RH) iy1 = RH-1;

        for (int py = iy0; py <= iy1; py++) {
            for (int px = ix0; px <= ix1; px++) {
                float qx = px + 0.5f, qy = py + 0.5f;
                float w0 = edge2(x0,y0,x1,y1,qx,qy);
                float w1 = edge2(x1,y1,x2,y2,qx,qy);
                float w2 = edge2(x2,y2,x0,y0,qx,qy);
                int in = (w0>=0 && w1>=0 && w2>=0) ||
                         (w0<=0 && w1<=0 && w2<=0);
                if (!in) continue;
                unsigned char* dst = &rgb[(py*RW+px)*3];
                dst[0]=(unsigned char)fminf(255,col[0]);
                dst[1]=(unsigned char)fminf(255,col[1]);
                dst[2]=(unsigned char)fminf(255,col[2]);
            }
        }
    }
}

/* count near-black pixels whose center lies in the analytic hole disk */
static void count_black_hole(const unsigned char* rgb, int* black, int* roi) {
    int bp = 0, rp = 0;
    for (int py = 0; py < RH; py++) {
        for (int px = 0; px < RW; px++) {
            float dx = ((px+0.5f) - RW*0.5f) / PX_PER_UNIT;
            float dy = -((py+0.5f) - RH*0.5f) / PX_PER_UNIT;
            if (dx*dx + dy*dy > ROI_FACTOR*ROI_FACTOR) continue;
            rp++;
            const unsigned char* p = &rgb[(py*RW+px)*3];
            float lum = 0.299f*p[0] + 0.587f*p[1] + 0.114f*p[2];
            if (lum < BLACK_LUM) bp++;
        }
    }
    *black = bp; *roi = rp;
}

static Metrics evaluate(const Mesh* m, const Loop* hole, int ci,
                        int fstart, int fcount, const Basis* b,
                        unsigned char* rgb) {
    Metrics M; memset(&M, 0, sizeof(M));
    M.dir_ratio = direction_ratio(m, fstart, fcount, b->n);
    M.dir_ok = (M.dir_ratio > 0.999f);
    M.outward_v = outward_patch_vertices(m, hole, ci, fstart, fcount, b->n);
    M.patch_v = hole->n + 1;
    M.nv_all = m->nv;
    outward_full_mesh(m, hole, ci, b->n, &M.outward_all,
                      &M.patch_out_allmesh);
    M.area_fill = fill_area(m, fstart, fcount);
    M.area_proj = loop_projected_area(m, hole);
    M.area_ratio = M.area_proj > 1e-12f ? M.area_fill / M.area_proj : 0.0f;
    M.area_ok = (M.area_ratio > 0.0f && M.area_ratio <= AREA_RATIO_MAX);
    render_mesh(m, fstart, fcount, b->u, b->v, b->n, rgb);
    count_black_hole(rgb, &M.black_px, &M.roi_px);
    return M;
}

/* ---------------- output artifacts ---------------- */

static void write_obj(const char* dir, const char* tag, const char* kind,
                      const Mesh* m) {
    char path[512];
    snprintf(path, sizeof(path), "%s/%s_%s.obj", dir, tag, kind);
    FILE* fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "cannot write %s\n", path); return; }
    fprintf(fp, "# holefill %s / %s ; outward CCW winding\n", tag, kind);
    for (int i = 0; i < m->nv; i++)
        fprintf(fp, "v %.6f %.6f %.6f\n", m->v[i].x, m->v[i].y, m->v[i].z);
    for (int i = 0; i < m->nf; i++)
        fprintf(fp, "f %d %d %d\n", m->f[i].a+1, m->f[i].b+1, m->f[i].c+1);
    fclose(fp);
}

static void write_ppm(const char* dir, const char* tag, const char* kind,
                      const unsigned char* rgb) {    char path[512];
    snprintf(path, sizeof(path), "%s/%s_%s.ppm", dir, tag, kind);
    FILE* fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "cannot write %s\n", path); return; }
    fprintf(fp, "P6\n%d %d\n255\n", RW, RH);
    fwrite(rgb, 1, (size_t)RW*RH*3, fp);
    fclose(fp);
}

static void write_json(const char* dir,
                       const Metrics good[3], const Metrics buggy[3],
                       const V3 norms[3]) {
    char path[512];
    snprintf(path, sizeof(path), "%s/holefill_metrics.json", dir);
    FILE* fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "cannot write %s\n", path); return; }
    fprintf(fp, "{\n");
    fprintf(fp, "  \"render\": {\"width\": %d, \"height\": %d, "
                "\"px_per_unit\": %.1f, \"roi_factor\": %.2f, "
                "\"black_lum_threshold\": %.0f},\n",
            RW, RH, PX_PER_UNIT, ROI_FACTOR, BLACK_LUM);
    fprintf(fp, "  \"area_ratio_max\": %.2f,\n", AREA_RATIO_MAX);
    fprintf(fp, "  \"fixed_reference_direction\": [0, 0, 1],\n");
    fprintf(fp, "  \"models\": [\n");
    const char* nm[3] = {"up","side","down"};
    for (int i = 0; i < 3; i++) {
        float dotref = vdot(norms[i], FIXED_UP);
        fprintf(fp,
            "    {\"name\": \"%s\", \"hole_normal\": [%.6f, %.6f, %.6f], "
            "\"ref_dot\": %.6f,\n"
            "      \"correct_ring_sign\": {\"dir_ratio\": %.3f, "
            "\"outward_vertices\": %d, \"patch_vertices\": %d, "
            "\"outward_allmesh\": %d, \"allmesh_vertices\": %d, "
            "\"patch_outward_allmesh\": %d, "
            "\"black_pixels\": %d, \"roi_pixels\": %d, "
            "\"area_ratio\": %.4f, \"area_ok\": %s, \"dir_ok\": %s},\n"
            "      \"buggy_fixed_up\": {\"dir_ratio\": %.3f, "
            "\"outward_vertices\": %d, \"outward_allmesh\": %d, "
            "\"patch_outward_allmesh\": %d, "
            "\"black_pixels\": %d, \"area_ratio\": %.4f}}%s\n",
            nm[i], norms[i].x, norms[i].y, norms[i].z, dotref,
            good[i].dir_ratio, good[i].outward_v, good[i].patch_v,
            good[i].outward_all, good[i].nv_all,
            good[i].patch_out_allmesh,
            good[i].black_px, good[i].roi_px,
            good[i].area_ratio, good[i].area_ok ? "true":"false",
            good[i].dir_ok ? "true":"false",
            buggy[i].dir_ratio, buggy[i].outward_v,
            buggy[i].outward_all, buggy[i].patch_out_allmesh,
            buggy[i].black_px, buggy[i].area_ratio,
            i < 2 ? "," : "");
    }
    fprintf(fp, "  ]\n}\n");
    fclose(fp);
}

/* ---------------- self checks / main ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    const char* outdir = argc > 1 ? argv[1] : NULL;
    if (outdir) mkdir(outdir, 0775);
    static unsigned char rgb[RW*RH*3];

    Metrics good[3], buggy[3];
    V3 norms[3];
    int all_struct_ok = 1;

    for (int k = 0; k < 3; k++) {
        Basis b = basis_of((HoleKind)k);
        norms[k] = b.n;
        char buf[96];

        Mesh m0 = build_flange((HoleKind)k);

        DEdge be[MAXBE];
        int ne = boundary_directed(&m0, be);
        Loop loops[MAXLOOPS];
        int nl = extract_loops(be, ne, loops);

        snprintf(buf, sizeof(buf), "%s_two_loops", b.name);
        int struct_ok = (ne == 12 && nl == 2);
        for (int i = 0; i < nl && struct_ok; i++)
            if (loops[i].n != 6) struct_ok = 0;
        check(buf, struct_ok);
        all_struct_ok &= struct_ok;

        Loop hole = select_hole_loop(&m0, loops, nl);
        snprintf(buf, sizeof(buf), "%s_hole_six", b.name);
        check(buf, hole.n == 6);

        float Aproj = loop_projected_area(&m0, &hole);
        snprintf(buf, sizeof(buf), "%s_proj_area", b.name);
        check(buf, fabsf(Aproj - 2.598076f) < 0.01f); /* hex r=1 */

        /* the decisive orientation-invariance property: the face-directed
         * ring signed normal is opposite truth for EVERY world orientation */
        V3 nL = vnorm(loop_newell(&m0, &hole));
        snprintf(buf, sizeof(buf), "%s_ring_sign_inward", b.name);
        check(buf, vdot(nL, b.n) < -0.999f);

        /* adjacent-edge cross sum carries the same sign as Newell */
        V3 nAdj = vnorm(loop_adjacent_cross(&m0, &hole));
        snprintf(buf, sizeof(buf), "%s_adjcross_matches_newell", b.name);
        check(buf, vdot(nAdj, nL) > 0.999f);

        /* every boundary owner tri is outward-wound vs truth */
        int owners_out = 1;
        for (int i = 0; i < ne; i++)
            if (vdot(face_raw_normal(&m0, be[i].owner), b.n) <= 0.0f)
                owners_out = 0;
        snprintf(buf, sizeof(buf), "%s_owners_outward", b.name);
        check(buf, owners_out);

        /* ---- correct ring-sign fill ---- */
        Mesh mc = m0;
        int fs, fc;
        int ci = fill_centroid_fan(&mc, &hole, ORIENT_RING_SIGN, &fs, &fc);
        good[k] = evaluate(&mc, &hole, ci, fs, fc, &b, rgb);

        snprintf(buf, sizeof(buf), "%s_dir_all_correct", b.name);
        check(buf, good[k].dir_ok);
        snprintf(buf, sizeof(buf), "%s_outward_seven", b.name);
        check(buf, good[k].outward_v == 7);
        snprintf(buf, sizeof(buf), "%s_allmesh_thirteen", b.name);
        check(buf, good[k].outward_all == 13 &&
                   good[k].patch_out_allmesh == 7);
        snprintf(buf, sizeof(buf), "%s_no_black_hole", b.name);
        check(buf, good[k].black_px == 0);
        snprintf(buf, sizeof(buf), "%s_area_gate_pass", b.name);
        check(buf, good[k].area_ok &&
                   fabsf(good[k].area_ratio - 1.0f) < 0.02f);
        if (outdir) {
            write_ppm(outdir, "correct", b.name, rgb);
            write_obj(outdir, "correct", b.name, &mc);
        }

        /* ---- buggy fixed-world-up fill on an independent copy ---- */
        Mesh mb = m0;
        int fb0, fbc;
        int cib = fill_centroid_fan(&mb, &hole, ORIENT_FIXED_UP, &fb0, &fbc);
        buggy[k] = evaluate(&mb, &hole, cib, fb0, fbc, &b, rgb);
        if (outdir) {
            write_ppm(outdir, "buggy", b.name, rgb);
            write_obj(outdir, "buggy", b.name, &mb);
        }
    }

    /* the fixed-reference diagnostic, asserted as the failure signature */
    check("buggy_up_still_ok",
          buggy[HOLE_UP].dir_ok && buggy[HOLE_UP].black_px == 0);
    check("buggy_down_all_black_normal_area",
          buggy[HOLE_DOWN].dir_ratio < 0.01f &&
          buggy[HOLE_DOWN].outward_v == 0 &&
          buggy[HOLE_DOWN].black_px == buggy[HOLE_DOWN].roi_px &&
          fabsf(buggy[HOLE_DOWN].area_ratio - 1.0f) < 0.02f);
    /* whole-mesh area weighting masks the flip at the ring: only the fan
     * center is inward (12/13); patch-scoped normals still say 0/7 */
    check("buggy_down_allmesh_masks_at_ring",
          buggy[HOLE_DOWN].outward_all == 12 &&
          buggy[HOLE_DOWN].patch_out_allmesh == 6);
    check("buggy_side_lost_reference_signal",
          fabsf(vdot(norms[HOLE_SIDE], FIXED_UP)) < 1e-6f &&
          buggy[HOLE_SIDE].dir_ratio < 0.01f);

    /* second criterion alone catches overshoot even when winding is right */
    {
        Basis b = basis_of(HOLE_UP);
        Mesh m = build_flange(HOLE_UP);
        DEdge be[MAXBE]; int ne = boundary_directed(&m, be);
        Loop loops[MAXLOOPS]; int nl = extract_loops(be, ne, loops);
        Loop hole = select_hole_loop(&m, loops, nl);
        int fs, fc;
        int ci = fill_balloon(&m, &hole, b.n, 2.0f, &fs, &fc);
        Metrics B = evaluate(&m, &hole, ci, fs, fc, &b, rgb);
        printf("balloon_area_ratio=%.4f\n", B.area_ratio);
        check("balloon_winding_still_right", B.dir_ok);
        check("balloon_area_gate_rejects",
              !B.area_ok && B.area_ratio > AREA_RATIO_MAX);
    }

    /* ---------------- report: three groups of four numbers ---------------- */
    printf("\n=== corrected fill (boundary-ring cross sign + area gate) ===\n");
    printf("%-5s %-26s %7s %8s %10s %9s %9s %12s\n",
           "model", "hole_normal", "n.ref", "dir_ok", "out_patch",
           "black/roi", "area_r", "out_allmesh");
    for (int k = 0; k < 3; k++) {
        printf("%-5s (%5.2f,%5.2f,%5.2f) %+7.3f %8.3f %5d/%-5d "
               "%9d/%-5d %9.4f %6d/%-5d\n",
               basis_of((HoleKind)k).name,
               norms[k].x, norms[k].y, norms[k].z,
               vdot(norms[k], FIXED_UP),
               good[k].dir_ratio, good[k].outward_v,
               good[k].patch_v, good[k].black_px, good[k].roi_px,
               good[k].area_ratio, good[k].outward_all, good[k].nv_all);
    }
    /* exactly 12 numbers: dir_ratio outward_v black_pixels area_ratio x3 */
    printf("TWELVE");
    for (int k = 0; k < 3; k++) {
        printf(" %.3f %d %d %.4f",
               good[k].dir_ratio, good[k].outward_v,
               good[k].black_px, good[k].area_ratio);
    }
    printf("\n");

    printf("\n=== buggy fill: fixed reference direction (0,0,1) ===\n");
    for (int k = 0; k < 3; k++) {
        printf("%-5s n=(%5.2f,%5.2f,%5.2f) n.ref=%+.0f  dir=%.3f out_patch=%d "
               "out_allmesh=%d/%d black=%d/%d area_r=%.4f\n",
               basis_of((HoleKind)k).name,
               norms[k].x, norms[k].y, norms[k].z,
               vdot(norms[k], FIXED_UP),
               buggy[k].dir_ratio, buggy[k].outward_v,
               buggy[k].outward_all, buggy[k].nv_all,
               buggy[k].black_px, buggy[k].roi_px, buggy[k].area_ratio);
    }
    printf("TWELVE_BUGGY");
    for (int k = 0; k < 3; k++) {
        printf(" %.3f %d %d %.4f",
               buggy[k].dir_ratio, buggy[k].outward_v,
               buggy[k].black_px, buggy[k].area_ratio);
    }
    printf("\n");

    printf("\nreference direction vs three hole normals:\n");
    printf("  up    n=( 0, 0,+1)  ref=(0,0,1)  dot=+1  angle=  0 deg\n");
    printf("  side  n=(+1, 0, 0)  ref=(0,0,1)  dot= 0  angle= 90 deg\n");
    printf("  down  n=( 0, 0,-1)  ref=(0,0,1)  dot=-1  angle=180 deg\n");

    if (outdir) {
        write_json(outdir, good, buggy, norms);
        printf("\nartifacts written under %s\n", outdir);
    }

    printf("\nall_ok=%d\n", failures == 0 && all_struct_ok ? 1 : 0);
    return (failures == 0 && all_struct_ok) ? 0 : 1;
}
