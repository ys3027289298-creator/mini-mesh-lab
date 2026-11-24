
/* lod: edge-collapse + position-weld simplification with UV-border locking.
 *
 * Background: LOD generation merges vertices by 3D position using a single
 * index stream (position and UV share an index, as in welded/). On a UV seam
 * two vertices coincide in 3D but carry different UVs; merging them re-indexes
 * every corner onto the surviving vertex, so its UV replaces the dropped one.
 * At a border, collapsing a boundary vertex into an interior one moves that
 * corner's UV inward; the adjacent UV triangle deforms outward across the
 * island margin and samples background texels. Up close the texture reads
 * background past the model boundary. Vertices on a UV seam/border are
 * therefore marked locked and may not be merged while simplifying.
 *
 * Lock rule (derived from island margin and texture size):
 *   margin_tex - outward clearance in texels: distance from the vertex UV to
 *                the texture border line, or to a non-incident UV boundary
 *                segment (background gap / facing edge of the other island),
 *                measured in UV space and scaled by the texture size;
 *   d_tex      - UV displacement one merge can cause, in texels. The physical
 *                collapse step is fixed by the 3D geometry (mean shortest
 *                incident 3D edge); multiply by the packed island's texel
 *                density (texels per 3D unit) so the gate stays comparable
 *                across margin tiers;
 *   weld jump  - UV distance to a coincident-position weld partner.
 * A boundary vertex is locked when its clearance cannot absorb the move:
 *   margin_tex <= d_tex,  or  margin_tex <= weld jump for any partner.
 * Smaller margins lock more vertices.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 128
#define MAXF 256
#define MAXE 512
#define TEXN 64
#define WELD_EPS 1e-4f

typedef struct { float x, y, z; } V3;
typedef struct { float u, v; } UV;
typedef struct { V3 p; UV uv; } Vert;
typedef struct { int a, b, c; } Tri;

typedef struct {
    Vert v[MAXV]; int nv;
    Tri f[MAXF]; int nf;
    int alive[MAXV];
} Mesh;

/* ---------- small helpers ---------- */

static float v3dist(V3 a, V3 b) {
    float dx = a.x-b.x, dy = a.y-b.y, dz = a.z-b.z;
    return sqrtf(dx*dx + dy*dy + dz*dz);
}
static float uvdist(UV a, UV b) {
    float du = a.u-b.u, dv = a.v-b.v;
    return sqrtf(du*du + dv*dv) * TEXN;        /* in texels */
}
static int same_v3(V3 a, V3 b) { return v3dist(a, b) < 1e-5f; }
static int same_uv(UV a, UV b) {
    return fabsf(a.u-b.u) < 1e-6f && fabsf(a.v-b.v) < 1e-6f;
}

static int mesh_add_vert(Mesh* m, V3 p, UV uv) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv].p = p; m->v[m->nv].uv = uv;
    m->alive[m->nv] = 1;
    return m->nv++;
}
static int mesh_add_tri(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a == b || b == c || a == c) return -1;
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* ---------- test geometry: gridded patch split into two UV islands ----------
 *
 * nx x ny quads in 3D; UVs split at column `split`. The split column gets
 * duplicated vertices (same 3D position, different UVs): the middle edge is an
 * internal UV seam, the outer ring a UV border. mo/gap in texels.
 */
static void build_grid(Mesh* m, int nx, int ny, int split,
                       float mo_tex, float gap_tex) {
    memset(m, 0, sizeof(*m));
    const float N = (float)TEXN;
    float lu0 = mo_tex/N, lu1 = (0.5f*N - 0.5f*gap_tex)/N;
    float ru0 = (0.5f*N + 0.5f*gap_tex)/N, ru1 = (N - mo_tex)/N;
    float vv0 = mo_tex/N, vv1 = (N - mo_tex)/N;

    int idL[32][32], idR[32][32];
    for (int ix = 0; ix <= nx; ix++) {
        for (int iy = 0; iy <= ny; iy++) {
            V3 p = { (float)ix, (float)iy, 0.0f };
            float vl = vv0 + (vv1-vv0) * ((float)iy / ny);
            float ul = lu0 + (lu1-lu0) * ((float)ix / split);
            float ur = ru0 + (ru1-ru0) * ((float)(ix-split) / (nx-split));
            if (ix < split)
                idL[ix][iy] = mesh_add_vert(m, p, (UV){ ul, vl });
            if (ix == split) {
                idL[ix][iy] = mesh_add_vert(m, p, (UV){ ul, vl });
                idR[ix][iy] = mesh_add_vert(m, p, (UV){ ur, vl }); /* seam duplicate */
            }
            if (ix > split)
                idR[ix][iy] = mesh_add_vert(m, p, (UV){ ur, vl });
        }
    }
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            int (*id)[32] = (ix < split) ? idL : idR;
            int a = id[ix][iy], b = id[ix+1][iy], c = id[ix][iy+1], d = id[ix+1][iy+1];
            mesh_add_tri(m, a, b, c);
            mesh_add_tri(m, b, d, c);
        }
    }
}

/* ---------- position-keyed edge records for seam detection ---------- */

typedef struct {
    int pa, pb;                 /* position ids */
    int va[4], vb[4];           /* corner vertex ids, up to 4 sides */
    UV  uva[4], uvb[4];
    int n;
} ERec;

static int position_groups(const Mesh* m, int posof[MAXV]) {
    int np = 0;
    for (int i = 0; i < m->nv; i++) {
        int g = -1;
        for (int k = 0; k < np && g < 0; k++)
            for (int j = 0; j < i; j++)
                if (posof[j] == k && same_v3(m->v[j].p, m->v[i].p)) { g = k; break; }
        posof[i] = (g < 0) ? np++ : g;
    }
    return np;
}

static void edge_records(const Mesh* m, const int posof[MAXV],
                         ERec rec[MAXE], int* nrec) {
    *nrec = 0;
    for (int t = 0; t < m->nf; t++) {
        int idx[3] = { m->f[t].a, m->f[t].b, m->f[t].c };
        for (int e = 0; e < 3; e++) {
            int va = idx[e], vb = idx[(e+1)%3];
            int pa = posof[va], pb = posof[vb];
            if (pa > pb) { int tmp = pa; pa = pb; pb = tmp; }
            int g = -1;
            for (int k = 0; k < *nrec; k++)
                if (rec[k].pa == pa && rec[k].pb == pb) { g = k; break; }
            if (g < 0) { g = (*nrec)++; rec[g].pa = pa; rec[g].pb = pb; rec[g].n = 0; }
            int s = rec[g].n++;
            rec[g].va[s] = va; rec[g].vb[s] = vb;
            rec[g].uva[s] = m->v[va].uv; rec[g].uvb[s] = m->v[vb].uv;
        }
    }
}

/* mesh border (one side) or internal seam (the two sides carry different UVs) */
static int edge_is_uv_boundary(const ERec* r) {
    if (r->n == 1) return 1;
    for (int s = 1; s < r->n; s++) {
        int match = (same_uv(r->uva[s], r->uva[0]) && same_uv(r->uvb[s], r->uvb[0])) ||
                    (same_uv(r->uva[s], r->uvb[0]) && same_uv(r->uvb[s], r->uva[0]));
        if (!match) return 1;
    }
    return 0;
}

/* ---------- geometric margin in texels ---------- */

static float point_seg_uv(UV p, UV a, UV b) {
    float ex = b.u-a.u, ey = b.v-a.v;
    float denom = ex*ex + ey*ey;
    float t = denom > 0.0f ? ((p.u-a.u)*ex + (p.v-a.v)*ey) / denom : 0.0f;
    if (t < 0.0f) t = 0.0f; else if (t > 1.0f) t = 1.0f;
    float cx = a.u + t*ex - p.u, cy = a.v + t*ey - p.v;
    return sqrtf(cx*cx + cy*cy) * TEXN;
}

/* outward clearance: texture border line, or nearest non-incident UV boundary */
static float margin_texels(const Mesh* m, int vi, const int posof[MAXV],
                           const ERec rec[MAXE], int nrec) {
    UV uv = m->v[vi].uv;
    float best = fminf(fminf(uv.u, 1.0f-uv.u), fminf(uv.v, 1.0f-uv.v)) * TEXN;
    int pi = posof[vi];
    for (int k = 0; k < nrec; k++) {
        if (!edge_is_uv_boundary(&rec[k])) continue;
        for (int s = 0; s < rec[k].n; s++) {
            int qa = posof[rec[k].va[s]], qb = posof[rec[k].vb[s]];
            if (qa == pi || qb == pi) continue;           /* incident edge */
            float d = point_seg_uv(uv, rec[k].uva[s], rec[k].uvb[s]);
            if (d < best) best = d;
        }
    }
    return best;
}

/* ---------- software rasterizer for bleed measurement ---------- */

static float edge_fn(float ax, float ay, float bx, float by, float px, float py) {
    return (px-ax)*(by-ay) - (py-ay)*(bx-ax);
}
static void raster_tri(const UV tri[3], int box[TEXN][TEXN], int value) {
    float ax = tri[0].u*TEXN, ay = tri[0].v*TEXN;
    float bx = tri[1].u*TEXN, by = tri[1].v*TEXN;
    float cx = tri[2].u*TEXN, cy = tri[2].v*TEXN;
    /* winding-independent inside test: sign follows the triangle's own area */
    float area = edge_fn(ax, ay, bx, by, cx, cy);
    if (fabsf(area) < 1e-9f) return;                  /* degenerate */
    float sign = area > 0.0f ? 1.0f : -1.0f;
    int x0 = (int)fminf(ax, fminf(bx, cx)); if (x0 < 0) x0 = 0;
    int y0 = (int)fminf(ay, fminf(by, cy)); if (y0 < 0) y0 = 0;
    int x1 = (int)ceilf(fmaxf(ax, fmaxf(bx, cx))); if (x1 > TEXN) x1 = TEXN;
    int y1 = (int)ceilf(fmaxf(ay, fmaxf(by, cy))); if (y1 > TEXN) y1 = TEXN;
    for (int py = y0; py < y1; py++)
        for (int px = x0; px < x1; px++) {
            float sx = px + 0.5f, sy = py + 0.5f;
            if (sign*edge_fn(bx,by,cx,cy,sx,sy) >= 0.0f &&
                sign*edge_fn(cx,cy,ax,ay,sx,sy) >= 0.0f &&
                sign*edge_fn(ax,ay,bx,by,sx,sy) >= 0.0f)
                box[py][px] = value;
        }
}

/* island labels painted from the original mesh: 0 background, 1/2 islands */
static void paint_islands(const Mesh* m, int label[TEXN][TEXN]) {
    memset(label, 0, sizeof(int)*TEXN*TEXN);
    for (int t = 0; t < m->nf; t++) {
        UV tri[3] = { m->v[m->f[t].a].uv, m->v[m->f[t].b].uv, m->v[m->f[t].c].uv };
        float cu = (tri[0].u + tri[1].u + tri[2].u) / 3.0f;
        raster_tri(tri, label, cu < 0.5f ? 1 : 2);
    }
}
static int bleed_pixels(const Mesh* m, const int label[TEXN][TEXN]) {
    int cover[TEXN][TEXN];
    memset(cover, 0, sizeof(cover));
    for (int t = 0; t < m->nf; t++) {
        UV tri[3] = { m->v[m->f[t].a].uv, m->v[m->f[t].b].uv, m->v[m->f[t].c].uv };
        raster_tri(tri, cover, 1);
    }
    int bleed = 0;
    for (int py = 0; py < TEXN; py++)
        for (int px = 0; px < TEXN; px++)
            if (cover[py][px] && label[py][px] == 0) bleed++;
    return bleed;
}

/* ---------- simplification: weld pairs + edge collapses ---------- */

typedef struct { int a, b; float cost; } Cand;

static int collect_candidates(const Mesh* m, Cand c[MAXE + MAXV]) {
    int n = 0;
    for (int i = 0; i < m->nv; i++) {
        if (!m->alive[i]) continue;
        for (int j = i+1; j < m->nv; j++) {
            if (!m->alive[j]) continue;
            if (v3dist(m->v[i].p, m->v[j].p) < WELD_EPS) {
                c[n].a = i; c[n].b = j; c[n].cost = 0.0f;   /* welds first */
                n++;
            }
        }
    }
    for (int t = 0; t < m->nf; t++) {
        int idx[3] = { m->f[t].a, m->f[t].b, m->f[t].c };
        for (int e = 0; e < 3; e++) {
            int a = idx[e], b = idx[(e+1)%3];
            if (a > b) { int tmp = a; a = b; b = tmp; }
            int seen = 0;
            for (int k = 0; k < n; k++)
                if (c[k].a == a && c[k].b == b) { seen = 1; break; }
            if (!seen) {
                c[n].a = a; c[n].b = b;
                c[n].cost = v3dist(m->v[a].p, m->v[b].p);
                n++;
            }
        }
    }
    return n;
}

/* Single-stream merge: corners are re-indexed onto `keep`, taking keep's UV
 * (the welded/ module does exactly this). A border vertex merged inward thus
 * drags its corner UVs and deforms the adjacent UV triangle across the margin.
 * Collapses move the kept position to the midpoint; welds leave it in place. */
static void merge_into(Mesh* m, int keep, int drop, int weld) {
    int nf = 0;
    for (int t = 0; t < m->nf; t++) {
        Tri r = m->f[t];
        if (r.a == drop) r.a = keep;
        if (r.b == drop) r.b = keep;
        if (r.c == drop) r.c = keep;
        if (r.a == r.b || r.b == r.c || r.a == r.c) continue;
        m->f[nf++] = r;
    }
    m->nf = nf;
    m->alive[drop] = 0;
    if (!weld) {
        m->v[keep].p.x = 0.5f*(m->v[keep].p.x + m->v[drop].p.x);
        m->v[keep].p.y = 0.5f*(m->v[keep].p.y + m->v[drop].p.y);
        m->v[keep].p.z = 0.5f*(m->v[keep].p.z + m->v[drop].p.z);
    }
}

static int simplify(Mesh* m, const int locked[MAXV], int target_faces) {
    while (m->nf > target_faces) {
        Cand c[MAXE + MAXV];
        int n = collect_candidates(m, c);
        int pick = -1;
        for (int k = 0; k < n; k++) {
            if (locked[c[k].a] || locked[c[k].b]) continue;
            if (pick < 0 || c[k].cost < c[pick].cost) pick = k;
        }
        if (pick < 0) break;
        merge_into(m, c[pick].a, c[pick].b, c[pick].cost == 0.0f);
    }
    return m->nf;
}

/* physical collapse step: mean per-vertex shortest incident 3D edge */
static float mean_min_edge3(const Mesh* m) {
    float sum = 0.0f; int used = 0;
    for (int i = 0; i < m->nv; i++) {
        if (!m->alive[i]) continue;
        float mn = 1e30f;
        for (int t = 0; t < m->nf; t++) {
            int idx[3] = { m->f[t].a, m->f[t].b, m->f[t].c };
            for (int e = 0; e < 3; e++)
                if (idx[e] == i) {
                    int j = idx[(e+1)%3];
                    float d = v3dist(m->v[i].p, m->v[j].p);
                    if (d < mn) mn = d;
                }
        }
        if (mn < 1e29f) { sum += mn; used++; }
    }
    return used ? sum / used : 0.0f;
}

/* largest UV jump from vi to any coincident-position weld partner */
static float max_weld_jump(const Mesh* m, int vi) {
    float best = 0.0f;
    for (int j = 0; j < m->nv; j++) {
        if (j == vi || !m->alive[j]) continue;
        if (v3dist(m->v[vi].p, m->v[j].p) < WELD_EPS) {
            float d = uvdist(m->v[vi].uv, m->v[j].uv);
            if (d > best) best = d;
        }
    }
    return best;
}

/* ---------- per-tier experiment ---------- */

typedef struct {
    const char* name;
    float mo, gap, threshold;
    int locked, seam_locked, border_locked;
    int bleed, faces, bleed_nolock;
} Tier;

static Tier run_tier(const char* name, float mo, float gap, int target_faces,
                     ERec seams_out[MAXE], float seam_margin[MAXE],
                     float seam_weld[MAXE], int* nseams) {
    Tier r; memset(&r, 0, sizeof(r));
    r.name = name; r.mo = mo; r.gap = gap;

    Mesh m;
    build_grid(&m, 8, 4, 4, mo, gap);

    int label[TEXN][TEXN];
    paint_islands(&m, label);

    int posof[MAXV];
    position_groups(&m, posof);
    ERec rec[MAXE]; int nrec;
    edge_records(&m, posof, rec, &nrec);

    int boundary[MAXV] = {0};
    int on_seam[MAXV] = {0};
    float margin[MAXV];
    for (int k = 0; k < nrec; k++) {
        if (!edge_is_uv_boundary(&rec[k])) continue;
        for (int s = 0; s < rec[k].n; s++) {
            boundary[rec[k].va[s]] = 1;
            boundary[rec[k].vb[s]] = 1;
            if (rec[k].n > 1) {                 /* internal UV seam */
                on_seam[rec[k].va[s]] = 1;
                on_seam[rec[k].vb[s]] = 1;
            }
        }
    }
    for (int i = 0; i < m.nv; i++)
        margin[i] = boundary[i] ? margin_texels(&m, i, posof, rec, nrec) : 1e30f;

    /* Reference texel density fixed at the widest-tier packing, so the gate
     * d_tex is identical across tiers and margin alone drives locking. */
    const float ref_mo = 8.0f, ref_gap = 16.0f;
    float texels_per_unit = (0.5f*TEXN - 0.5f*ref_gap - ref_mo) / 4.0f;
    float step = mean_min_edge3(&m) * texels_per_unit;
    r.threshold = step;

    int locked[MAXV] = {0};
    for (int i = 0; i < m.nv; i++) {
        if (!boundary[i]) continue;
        if (on_seam[i]) {
            /* A seam vertex has a coincident-position partner on the other
             * island. A single weld re-indexes across the whole inter-island
             * gap (weld jump), which no outward margin can absorb, so seam
             * vertices are locked unconditionally. */
            locked[i] = 1;
        } else {
            /* outer border: lock when the clearance cannot absorb one step */
            locked[i] = (margin[i] <= step);
        }
        if (locked[i]) {
            if (on_seam[i]) r.seam_locked++; else r.border_locked++;
        }
    }
    for (int i = 0; i < m.nv; i++) if (locked[i]) r.locked++;

    /* Seam diagnostic snapshot taken BEFORE any merge: margins and weld jumps
     * must be measured on the intact mesh, while every seam partner is alive. */
    *nseams = 0;
    for (int k = 0; k < nrec; k++) {
        if (!edge_is_uv_boundary(&rec[k]) || rec[k].n == 1) continue;  /* seam only */
        if (*nseams >= MAXE) break;
        seams_out[*nseams] = rec[k];
        float sm = 1e30f, wj = 0.0f;
        for (int s = 0; s < rec[k].n; s++) {
            int va = rec[k].va[s];
            if (margin[va] < sm) sm = margin[va];
            float w = max_weld_jump(&m, va);
            if (w > wj) wj = w;
        }
        seam_margin[*nseams] = sm;
        seam_weld[*nseams] = wj;
        (*nseams)++;
    }

    /* snapshot locked geometry to prove it never moves */
    V3 lock_p[MAXV]; UV lock_uv[MAXV];
    for (int i = 0; i < m.nv; i++) { lock_p[i] = m.v[i].p; lock_uv[i] = m.v[i].uv; }

    /* control: same simplification with every lock disabled (the old bug) */
    Mesh nolock = m;
    simplify(&nolock, (int[MAXV]){0}, target_faces);
    r.bleed_nolock = bleed_pixels(&nolock, label);

    r.faces = simplify(&m, locked, target_faces);
    r.bleed = bleed_pixels(&m, label);

    /* verify every locked vertex survived unmoved */
    for (int i = 0; i < m.nv; i++) {
        if (!locked[i]) continue;
        if (!(m.alive[i] && same_v3(lock_p[i], m.v[i].p) && same_uv(lock_uv[i], m.v[i].uv)))
            r.locked = -1;
    }
    return r;
}

/* ---------- self tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const int target = 32;     /* halve the 64 input faces */

    Mesh probe;
    build_grid(&probe, 8, 4, 4, 8.0f, 16.0f);
    printf("probe: nv=%d nf=%d\n", probe.nv, probe.nf);
    check("grid_verts", probe.nv == 45 + 5);     /* 45 positions + 5 seam duplicates */
    check("grid_faces", probe.nf == 64);

    int label[TEXN][TEXN];
    paint_islands(&probe, label);
    check("unpainted_no_bleed", bleed_pixels(&probe, label) == 0);

    int posof[MAXV];
    position_groups(&probe, posof);
    ERec rec[MAXE]; int nrec;
    edge_records(&probe, posof, rec, &nrec);
    int nb = 0, nseam = 0;
    for (int k = 0; k < nrec; k++)
        if (edge_is_uv_boundary(&rec[k])) { nb++; if (rec[k].n > 1) nseam++; }
    printf("probe: position edges=%d uv boundary=%d seams=%d\n", nrec, nb, nseam);
    check("uv_boundary_edges", nb == 2*(8+4) + 4);   /* outer ring + 4 seam edges */
    check("uv_seams", nseam == 4);

    ERec seams[MAXE]; float seam_margins[MAXE], seam_welds[MAXE]; int nseams;
    Tier wide   = run_tier("wide",   8.0f, 16.0f, target, seams, seam_margins, seam_welds, &nseams);
    Tier medium = run_tier("medium", 4.0f,  6.0f, target, seams, seam_margins, seam_welds, &nseams);
    Tier narrow = run_tier("narrow", 1.0f,  2.0f, target, seams, seam_margins, seam_welds, &nseams);

    check("locks_held", wide.locked >= 0 && medium.locked >= 0 && narrow.locked >= 0);
    check("lock_monotonic",
          narrow.locked >= medium.locked && medium.locked >= wide.locked);
    check("bug_reproduced_nolock",
          wide.bleed_nolock > 0 && medium.bleed_nolock > 0 && narrow.bleed_nolock > 0);
    check("locked_runs_clean",
          wide.bleed == 0 && medium.bleed == 0 && narrow.bleed == 0);
    check("faces_reduced",
          wide.faces <= 64 && medium.faces <= 64 && narrow.faces <= 64);

    printf("\n== LOD experiment: texture %dx%d, target %d faces (input 64) ==\n",
           TEXN, TEXN, target);
    printf("%-7s %8s %8s %10s %9s %8s %8s %14s\n",
           "tier", "mo_tex", "gap_tex", "thresh_tex", "locked(s+b)", "bleed",
           "faces", "bleed_no_lock");
    const Tier tiers[3] = { wide, medium, narrow };
    for (int i = 0; i < 3; i++)
        printf("%-7s %8.1f %8.1f %10.2f %5d(%d+%d) %8d %8d %14d\n",
               tiers[i].name, tiers[i].mo, tiers[i].gap, tiers[i].threshold,
               tiers[i].locked, tiers[i].seam_locked, tiers[i].border_locked,
               tiers[i].bleed, tiers[i].faces,
               tiers[i].bleed_nolock);

    printf("\nnine numbers (locked, bleed, faces):\n");
    for (int i = 0; i < 3; i++)
        printf("%s: %d %d %d\n", tiers[i].name,
               tiers[i].locked, tiers[i].bleed, tiers[i].faces);

    if (narrow.bleed > 0) {
        /* lock criterion missed a UV seam: put seam margin and weld jump next
         * to the gate so the missed case is visible */
        printf("\nDIAGNOSTIC narrow tier: residual bleed=%d\n", narrow.bleed);
        printf("lock gate: d_tex=%.2f texels; per-seam clearance vs gate:\n", narrow.threshold);
        printf("%-8s %12s %12s %12s %10s\n",
               "seam", "margin_tex", "weldjump_tex", "thresh_tex", "verdict");
        for (int k = 0; k < nseams; k++)
            printf("p%02d-p%02d %12.2f %12.2f %12.2f %10s\n",
                   seams[k].pa, seams[k].pb, seam_margins[k], seam_welds[k],
                   narrow.threshold,
                   seam_margins[k] <= narrow.threshold ? "locked" : "MISSED");
    } else {
        printf("\nnarrow tier bleed=0: every UV seam was covered by the lock "
               "criterion; diagnostic not triggered.\n");
    }

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
