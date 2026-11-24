/* seam: tolerance-driven normal merging across UV seams. */
#include "seam.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- folded-strip generator constants (unscaled units; everything x scale) ---- */
#define PANELS        SEAM_PANELS
#define COL_PER_PANEL SEAM_COL_PER_PANEL
#define ROW           SEAM_ROW
#define CELL          SEAM_CELL
#define F_SEAM        SEAM_F_SEAM
#define F_CREASE      SEAM_F_CREASE
#define BUMP          0.10f   /* smooth panel curvature amplitude */
#define SLOPE_ABS     1.00f   /* alternating panel slopes -> ~46 deg fold */

static V3 vsub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 vadd(V3 p, V3 q) { V3 r = { p.x+q.x, p.y+q.y, p.z+q.z }; return r; }
static V3 vcross(V3 p, V3 q) {
    V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r;
}
static float vdot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float vlen(V3 p) { return sqrtf(vdot(p, p)); }
static V3 vnorm(V3 p) {
    float l = vlen(p); V3 r = {0,0,0};
    if (l > 1e-12f) { r.x=p.x/l; r.y=p.y/l; r.z=p.z/l; }
    return r;
}

void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

int mesh_add_vertex(Mesh* m, V3 p, UV uv, int seam_flag) {
    if (m->nv >= MAXV) return -1;
    m->pos[m->nv] = p;
    m->uv[m->nv] = uv;
    m->seam[m->nv] = seam_flag;
    return m->nv++;
}

int mesh_add_tri(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a < 0 || b < 0 || c < 0 || a >= m->nv || b >= m->nv || c >= m->nv) return -1;
    if (a == b || b == c || a == c) return -1;
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* panel surface: height above the fold line at local coordinate q in [0,1] */
static float panel_height(int k, float q) {
    float h = 0.0f;
    for (int j = 0; j < k; j++) h += (j % 2 == 0 ? SLOPE_ABS : -SLOPE_ABS);
    float s = (k % 2 == 0 ? SLOPE_ABS : -SLOPE_ABS);
    return h + s * q + BUMP * sinf((float)M_PI * q);
}
static float panel_dydx(int k, float q) {
    float s = (k % 2 == 0 ? SLOPE_ABS : -SLOPE_ABS);
    return s + BUMP * (float)M_PI * cosf((float)M_PI * q);
}
/* analytic smooth normal of panel k at local q */
static V3 panel_normal(int k, float q) {
    V3 n = { -panel_dydx(k, q), 1.0f, 0.0f };
    return vnorm(n);
}

/* choose nseam interior grid lines (global interior line index 0..PANELS*(COL_PER_PANEL-1)-1) */
static void pick_seam_lines(int nseam, int* lines) {
    int total = PANELS * (COL_PER_PANEL - 1);
    for (int i = 0; i < nseam; i++) {
        lines[i] = (int)((i + 0.5f) * total / nseam);
        if (lines[i] >= total) lines[i] = total - 1;
    }
}

void mesh_build_folded_strip(Mesh* m, int nseam, float scale) {
    mesh_init(m);
    int seam_lines[COL_PER_PANEL * PANELS];
    pick_seam_lines(nseam, seam_lines);
    int is_seam[COL_PER_PANEL * PANELS];
    memset(is_seam, 0, sizeof(is_seam));
    for (int i = 0; i < nseam; i++) is_seam[seam_lines[i]] = 1;

    float crease_dx = F_CREASE * CELL * scale;
    float seam_dx   = F_SEAM   * CELL * scale;

    for (int k = 0; k < PANELS; k++) {
        float xoff = k * (scale + crease_dx);
        /* split the panel's column range into UV-island segments at the
         * seam lines that fall inside this panel.  Each segment allocates
         * its own vertex records, so adjacent segments share a world
         * position at the seam but never an index -> separate channels. */
        int starts[COL_PER_PANEL + 1];
        int nseg = 0;
        starts[nseg++] = 0;
        for (int c = 1; c < COL_PER_PANEL; c++) {
            int global_line = k * (COL_PER_PANEL - 1) + (c - 1);
            if (is_seam[global_line]) starts[nseg++] = c;
        }

        for (int sg = 0; sg < nseg; sg++) {
            int c0 = starts[sg];
            int c1 = (sg + 1 < nseg) ? starts[sg + 1] : COL_PER_PANEL;
            int base = m->nv;
            /* vertex grid for this segment, columns major:
             * record id = base + (c-c0)*(ROW+1) + r */
            #define SEG_ID(cc, rr) (base + ((cc) - c0) * (ROW + 1) + (rr))
            for (int c = c0; c <= c1; c++) {
                /* a record is seam-flagged only when this column is an
                 * interior UV split inside the panel; panel rims (c==0,
                 * c==COL_PER_PANEL) are crease boundaries, never seams */
                int interior = (c > 0 && c < COL_PER_PANEL);
                int at_seam = 0;
                if (interior) {
                    int gl = k * (COL_PER_PANEL - 1) + (c - 1);
                    at_seam = is_seam[gl];
                }
                for (int r = 0; r <= ROW; r++) {
                    float q = (float)c * CELL;
                    /* the segment starting after a seam shifts its first
                     * column by seam_dx; the segment ending at a seam
                     * keeps the nominal position */
                    float shift = (at_seam && c == c0) ? seam_dx : 0.0f;
                    float x = xoff + c * CELL * scale + shift;
                    float y = panel_height(k, q) * scale
                              + (shift > 0.0f ? panel_dydx(k, q) * seam_dx : 0.0f);
                    float z = (float)r / ROW * scale;
                    V3 p = { x, y, z };
                    UV uv = { (float)c / COL_PER_PANEL + shift / scale,
                              (float)r / ROW };
                    int id = mesh_add_vertex(m, p, uv, at_seam);
                    m->truth[id] = panel_normal(k, q);
                }
            }
            for (int c = c0; c < c1; c++) {
                for (int r = 0; r < ROW; r++) {
                    mesh_add_tri(m, SEG_ID(c, r), SEG_ID(c + 1, r + 1), SEG_ID(c + 1, r));
                    mesh_add_tri(m, SEG_ID(c, r), SEG_ID(c, r + 1), SEG_ID(c + 1, r + 1));
                }
            }
            #undef SEG_ID
        }
    }
    mesh_normals(m);
}

float mesh_characteristic_length(const Mesh* m) {
    if (m->nf == 0) return -1.0f;
    int cap = m->nf * 3;
    float* e = malloc(sizeof(float) * (size_t)cap);
    if (!e) return -1.0f;
    int n = 0;
    for (int i = 0; i < m->nf; i++) {
        int idx[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
        for (int j = 0; j < 3; j++) {
            e[n++] = vlen(vsub(m->pos[idx[j]], m->pos[idx[(j + 1) % 3]]));
        }
    }
    /* insertion sort median (n is a few thousand) */
    for (int i = 1; i < n; i++) {
        float key = e[i]; int j = i - 1;
        while (j >= 0 && e[j] > key) { e[j + 1] = e[j]; j--; }
        e[j + 1] = key;
    }
    float med = e[n / 2];
    free(e);
    return med;
}

void mesh_normals(Mesh* m) {
    for (int i = 0; i < m->nv; i++) { m->nrm[i].x = m->nrm[i].y = m->nrm[i].z = 0.0f; }
    for (int i = 0; i < m->nf; i++) {
        V3 a = m->pos[m->f[i].a], b = m->pos[m->f[i].b], c = m->pos[m->f[i].c];
        V3 n = vcross(vsub(b, a), vsub(c, a));
        int idx[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
        for (int j = 0; j < 3; j++) {
            m->nrm[idx[j]] = vadd(m->nrm[idx[j]], n);
        }
    }
    for (int i = 0; i < m->nv; i++) m->nrm[i] = vnorm(m->nrm[i]);
}

void mesh_ground_truth_normals(const Mesh* m, V3* out) {
    for (int i = 0; i < m->nv; i++) out[i] = m->truth[i];
}

void mesh_find_pairs(const Mesh* m, PairSets* ps) {
    memset(ps, 0, sizeof(*ps));
    /* recover the build scale from the x extent:
     * extent = PANELS*scale + (PANELS-1)*F_CREASE*CELL*scale */
    float ext = 0.0f;
    for (int i = 0; i < m->nv; i++)
        if (m->pos[i].x > ext) ext = m->pos[i].x;
    float scale = ext / (PANELS + (PANELS - 1) * F_CREASE * CELL);
    float r = 0.1f * scale;   /* only true twins are closer than 0.1 cell */
    for (int i = 0; i < m->nv; i++) {
        for (int j = i + 1; j < m->nv; j++) {
            V3 d3 = vsub(m->pos[i], m->pos[j]);
            float d = vlen(d3);
            if (d > r || d <= 0.0f) continue;
            if (fabsf(m->pos[i].z - m->pos[j].z) > 1e-6f * scale) continue;
            Pair pr = { i, j, d };
            if (m->seam[i] && m->seam[j]) {
                if (ps->nseam < MAXPAIRS) ps->seam[ps->nseam++] = pr;
            } else if (!m->seam[i] && !m->seam[j]) {
                if (ps->ncrease < MAXPAIRS) ps->crease[ps->ncrease++] = pr;
            }
        }
    }
}

/* ---- union-find ---- */
void uf_init(UF* uf, int n) { uf->n = n; for (int i = 0; i < n; i++) uf->parent[i] = i; }
int uf_find(UF* uf, int x) {
    while (uf->parent[x] != x) { uf->parent[x] = uf->parent[uf->parent[x]]; x = uf->parent[x]; }
    return x;
}
void uf_union(UF* uf, int a, int b) {
    a = uf_find(uf, a); b = uf_find(uf, b);
    if (a != b) {
        if (a < b) uf->parent[b] = a;
        else       uf->parent[a] = b;
    }
}

/* ---- spatial hash grid for neighbourhood welds ---- */
typedef struct { int* head; int* next; int cap; float inv; V3 origin; } Grid;

static void grid_build(Grid* g, const Mesh* m, float eps) {
    g->cap = 8;
    while (g->cap < m->nv * 2) g->cap <<= 1;
    g->head = malloc(sizeof(int) * (size_t)g->cap);
    g->next = malloc(sizeof(int) * (size_t)m->nv);
    for (int i = 0; i < g->cap; i++) g->head[i] = -1;
    g->inv = 1.0f / eps;
    g->origin = m->pos[0];
    for (int i = 0; i < m->nv; i++) {
        V3 d = vsub(m->pos[i], g->origin);
        int cx = (int)floorf(d.x * g->inv), cy = (int)floorf(d.y * g->inv),
            cz = (int)floorf(d.z * g->inv);
        unsigned h = (unsigned)(cx * 73856093 ^ cy * 19349663 ^ cz * 83492791) & (unsigned)(g->cap - 1);
        g->next[i] = g->head[h];
        g->head[h] = i;
    }
}
static void grid_free(Grid* g) { free(g->head); free(g->next); }

static void grid_weld(const Mesh* m, float eps, float angle_limit_deg, int guarded, UF* uf) {
    Grid g; grid_build(&g, m, eps);
    float cos_limit = cosf(angle_limit_deg * (float)M_PI / 180.0f);
    for (int i = 0; i < m->nv; i++) {
        V3 d = vsub(m->pos[i], g.origin);
        int cx = (int)floorf(d.x * g.inv), cy = (int)floorf(d.y * g.inv),
            cz = (int)floorf(d.z * g.inv);
        for (int ox = -1; ox <= 1; ox++)
        for (int oy = -1; oy <= 1; oy++)
        for (int oz = -1; oz <= 1; oz++) {
            unsigned h = (unsigned)((cx+ox) * 73856093 ^ (cy+oy) * 19349663 ^ (cz+oz) * 83492791)
                         & (unsigned)(g.cap - 1);
            for (int j = g.head[h]; j >= 0; j = g.next[j]) {
                if (j >= i) continue;
                if (vlen(vsub(m->pos[i], m->pos[j])) > eps) continue;
                if (guarded && vdot(m->nrm[i], m->nrm[j]) < cos_limit) continue;
                uf_union(uf, i, j);
            }
        }
    }
    grid_free(&g);
}

void weld_position(const Mesh* m, float eps, UF* uf) {
    uf_init(uf, m->nv);
    grid_weld(m, eps, 0.0f, 0, uf);
}
void weld_guarded(const Mesh* m, float eps, float angle_limit_deg, UF* uf) {
    uf_init(uf, m->nv);
    grid_weld(m, eps, angle_limit_deg, 1, uf);
}

void share_group_normals(const Mesh* m, const UF* uf, V3* out) {
    V3 acc[MAXV];
    for (int i = 0; i < m->nv; i++) { acc[i].x = acc[i].y = acc[i].z = 0.0f; }
    /* m->nrm are unit per-record normals; accumulate unit normals per group */
    for (int i = 0; i < m->nv; i++) {
        int r = uf_find((UF*)uf, i);
        acc[r] = vadd(acc[r], m->nrm[i]);
    }
    for (int i = 0; i < m->nv; i++) acc[i] = vnorm(acc[i]);
    for (int i = 0; i < m->nv; i++) out[i] = acc[uf_find((UF*)uf, i)];
}

int count_uv_channels(const Mesh* m) {
    UF uf; uf_init(&uf, m->nv);
    for (int i = 0; i < m->nf; i++) {
        uf_union(&uf, m->f[i].a, m->f[i].b);
        uf_union(&uf, m->f[i].b, m->f[i].c);
        uf_union(&uf, m->f[i].a, m->f[i].c);
    }
    int ch = 0;
    for (int i = 0; i < m->nv; i++) if (uf_find(&uf, i) == i) ch++;
    return ch;
}

int count_open_normal_channels(const Mesh* m, const UF* uf, const PairSets* ps) {
    (void)m;
    int open = 0;
    for (int i = 0; i < ps->nseam; i++)
        if (uf_find((UF*)uf, ps->seam[i].a) != uf_find((UF*)uf, ps->seam[i].b)) open++;
    return open;
}

/* seam pairs whose records are in different position-weld groups */
int count_unmerged_seams(const PairSets* ps, const UF* uf) {
    int open = 0;
    for (int i = 0; i < ps->nseam; i++)
        if (uf_find((UF*)uf, ps->seam[i].a) != uf_find((UF*)uf, ps->seam[i].b)) open++;
    return open;
}

int count_mismerged_creases(const PairSets* ps, const UF* uf) {
    int bad = 0;
    for (int i = 0; i < ps->ncrease; i++)
        if (uf_find((UF*)uf, ps->crease[i].a) == uf_find((UF*)uf, ps->crease[i].b)) bad++;
    return bad;
}

static int normal_equal(V3 a, V3 b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

/* seam pairs still carrying distinct normals after the normal pass */
int count_open_normal_seams(const PairSets* ps, const V3* nrm) {
    int open = 0;
    for (int i = 0; i < ps->nseam; i++)
        if (!normal_equal(nrm[ps->seam[i].a], nrm[ps->seam[i].b])) open++;
    return open;
}

/* distinct normal values per weld group, summed over groups */
int count_normal_values(const Mesh* m, const UF* uf, const V3* nrm) {
    int total = 0;
    for (int g = 0; g < m->nv; g++) {
        if (uf_find((UF*)uf, g) != g) continue;
        int distinct = 0;
        for (int i = 0; i < m->nv; i++) {
            if (uf_find((UF*)uf, i) != g) continue;
            int seen = 0;
            for (int j = 0; j < i; j++) {
                if (uf_find((UF*)uf, j) == g && normal_equal(nrm[i], nrm[j])) { seen = 1; break; }
            }
            if (!seen) distinct++;
        }
        total += distinct;
    }
    return total;
}

float mean_normal_error_deg(const Mesh* m, const V3* out) {
    double sum = 0.0;
    for (int i = 0; i < m->nv; i++) {
        float c = vdot(vnorm(out[i]), m->truth[i]);
        if (c > 1.0f) c = 1.0f;
        if (c < -1.0f) c = -1.0f;
        sum += acosf(c) * 180.0f / (float)M_PI;
    }
    return (float)(sum / m->nv);
}

/* ---- orthographic software rasteriser ---- */
#define CAM_ANGLE 25.0f
#define AMBIENT   0.35f

static V3 light_dir(void) {
    V3 L = { 0.40f, 0.80f, 0.45f };
    return vnorm(L);
}

typedef struct {
    float sx, sy, depth, shade;
    int seam;
} VScreen;

static VScreen project(V3 p, V3 n, int seam, float cam_c, float cam_s,
                       float s2w, V3 center) {
    VScreen r;
    r.sx = (p.x - center.x) * s2w + RENDER_W * 0.5f;
    r.sy = RENDER_H * 0.5f - ((p.z - center.z) * cam_c + (p.y - center.y) * cam_s) * s2w;
    r.depth = (p.y - center.y) * cam_c - (p.z - center.z) * cam_s;
    float lam = vdot(n, light_dir());
    r.shade = (lam < 0.0f ? 0.0f : lam) * (1.0f - AMBIENT) + AMBIENT;
    r.seam = seam;
    return r;
}

/* paint a thick screen line (Bresenham plus a 1px dilation) into mask */
static void paint_line(int* mask, float x0, float y0, float x1, float y1) {
    int ix0 = (int)lroundf(x0), iy0 = (int)lroundf(y0);
    int ix1 = (int)lroundf(x1), iy1 = (int)lroundf(y1);
    int dx = abs(ix1 - ix0), dy = -abs(iy1 - iy0);
    int sx = ix0 < ix1 ? 1 : -1, sy = iy0 < iy1 ? 1 : -1;
    int err = dx + dy;
    for (;;) {
        for (int oy = -1; oy <= 1; oy++)
            for (int ox = -1; ox <= 1; ox++) {
                int px = ix0 + ox, py = iy0 + oy;
                if (px >= 0 && px < RENDER_W && py >= 0 && py < RENDER_H)
                    mask[py * RENDER_W + px] = 1;
            }
        if (ix0 == ix1 && iy0 == iy1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; ix0 += sx; }
        if (e2 <= dx) { err += dx; iy0 += sy; }
    }
}

/* mask the exact screen projection of the seam lines.  Seam twins share
 * an (x,y) position and differ only in z-row, so group seam records into
 * rows of coincident (x,y), then join consecutive rows. */
static void build_seam_mask(const Mesh* m, const V3* nrm, int* mask,
                            float cam_c, float cam_s, float s2w, V3 center) {
    VScreen* sp = malloc(sizeof(VScreen) * (size_t)m->nv);
    if (!sp) return;
    for (int i = 0; i < m->nv; i++)
        sp[i] = project(m->pos[i], nrm[i], m->seam[i], cam_c, cam_s, s2w, center);
    float r2 = (0.5f * CELL * s2w) * (0.5f * CELL * s2w);
    for (int i = 0; i < m->nv; i++) {
        if (!m->seam[i]) continue;
        for (int j = i + 1; j < m->nv; j++) {
            if (!m->seam[j]) continue;
            float dx = sp[i].sx - sp[j].sx, dy = sp[i].sy - sp[j].sy;
            if (dx * dx + dy * dy > r2) continue;
            /* coincident column: only connect if the rows are adjacent in z */
            float dz = fabsf(m->pos[i].z - m->pos[j].z);
            float rowstep = 1.0f / ROW;
            if (dz < 1.5f * rowstep)
                paint_line(mask, sp[i].sx, sp[i].sy, sp[j].sx, sp[j].sy);
        }
    }
    free(sp);
}

static void rasterize(const Mesh* m, const V3* nrm, float* shade, int* mask, float* depth) {
    float cam_c = cosf(CAM_ANGLE * (float)M_PI / 180.0f);
    float cam_s = sinf(CAM_ANGLE * (float)M_PI / 180.0f);
    V3 mn = m->pos[0], mx = m->pos[0];
    for (int i = 1; i < m->nv; i++) {
        mn.x = fminf(mn.x, m->pos[i].x); mn.y = fminf(mn.y, m->pos[i].y); mn.z = fminf(mn.z, m->pos[i].z);
        mx.x = fmaxf(mx.x, m->pos[i].x); mx.y = fmaxf(mx.y, m->pos[i].y); mx.z = fmaxf(mx.z, m->pos[i].z);
    }
    V3 center = { (mn.x + mx.x) * 0.5f, (mn.y + mx.y) * 0.5f, (mn.z + mx.z) * 0.5f };
    float span_x = mx.x - mn.x;
    float span_s = (mx.z - mn.z) * cam_c + (mx.y - mn.y) * cam_s;
    float s2w = fminf(RENDER_W * 0.92f / span_x, RENDER_H * 0.92f / span_s);
    int N = RENDER_W * RENDER_H;
    for (int i = 0; i < N; i++) { shade[i] = 0.0f; mask[i] = 0; depth[i] = 1e30f; }

    for (int t = 0; t < m->nf; t++) {
        int ia = m->f[t].a, ib = m->f[t].b, ic = m->f[t].c;
        VScreen va = project(m->pos[ia], nrm[ia], m->seam[ia], cam_c, cam_s, s2w, center);
        VScreen vb = project(m->pos[ib], nrm[ib], m->seam[ib], cam_c, cam_s, s2w, center);
        VScreen vc = project(m->pos[ic], nrm[ic], m->seam[ic], cam_c, cam_s, s2w, center);
        float area = (vb.sx - va.sx) * (vc.sy - va.sy) - (vb.sy - va.sy) * (vc.sx - va.sx);
        if (fabsf(area) < 1e-9f) continue;
        int x0 = (int)fminf(va.sx, fminf(vb.sx, vc.sx));
        int x1 = (int)fmaxf(va.sx, fmaxf(vb.sx, vc.sx));
        int y0 = (int)fminf(va.sy, fminf(vb.sy, vc.sy));
        int y1 = (int)fmaxf(va.sy, fmaxf(vb.sy, vc.sy));
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 >= RENDER_W) x1 = RENDER_W - 1;
        if (y1 >= RENDER_H) y1 = RENDER_H - 1;
        for (int py = y0; py <= y1; py++) {
            for (int px = x0; px <= x1; px++) {
                float px_ = px + 0.5f, py_ = py + 0.5f;
                float w0 = ((vb.sx - px_) * (vc.sy - py_) - (vb.sy - py_) * (vc.sx - px_)) / area;
                float w1 = ((vc.sx - px_) * (va.sy - py_) - (vc.sy - py_) * (va.sx - px_)) / area;
                float w2 = 1.0f - w0 - w1;
                if (w0 < 0 || w1 < 0 || w2 < 0) continue;
                int idx = py * RENDER_W + px;
                float d = w0 * va.depth + w1 * vb.depth + w2 * vc.depth;
                if (d >= depth[idx]) continue;
                depth[idx] = d;
                shade[idx] = w0 * va.shade + w1 * vb.shade + w2 * vc.shade;
            }
        }
    }
    /* seam lines are a screen-space overlay, independent of the normal field */
    build_seam_mask(m, nrm, mask, cam_c, cam_s, s2w, center);
}

static void write_ppm(const char* path, const float* shade) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { perror(path); return; }
    fprintf(fp, "P6\n%d %d\n255\n", RENDER_W, RENDER_H);
    int N = RENDER_W * RENDER_H;
    unsigned char* row = malloc((size_t)N * 3);
    for (int i = 0; i < N; i++) {
        float s = shade[i];
        int g = (int)(s * 255.0f + 0.5f);
        if (g < 0) g = 0;
        if (g > 255) g = 255;
        row[i * 3] = (unsigned char)g;
        row[i * 3 + 1] = (unsigned char)g;
        row[i * 3 + 2] = (unsigned char)g;
    }
    fwrite(row, 1, (size_t)N * 3, fp);
    free(row);
    fclose(fp);
}

int render_bright_pixels(const Mesh* m, const V3* nrm, float shade_tol, const char* ppm_path) {
    int N = RENDER_W * RENDER_H;
    float *test = malloc(sizeof(float) * N), *truth = malloc(sizeof(float) * N),
          *dep1 = malloc(sizeof(float) * N), *dep2 = malloc(sizeof(float) * N);
    int* mask = malloc(sizeof(int) * N);
    if (!test || !truth || !dep1 || !dep2 || !mask) {
        free(test); free(truth); free(dep1); free(dep2); free(mask);
        return -1;
    }
    rasterize(m, nrm, test, mask, dep1);
    V3* tn = malloc(sizeof(V3) * (size_t)m->nv);
    mesh_ground_truth_normals(m, tn);
    rasterize(m, tn, truth, mask, dep2);
    int bright = 0;
    for (int i = 0; i < N; i++) {
        if (mask[i] && fabsf(test[i] - truth[i]) > shade_tol) bright++;
    }
    if (ppm_path) write_ppm(ppm_path, test);
    free(test); free(truth); free(dep1); free(dep2); free(mask); free(tn);
    return bright;
}
