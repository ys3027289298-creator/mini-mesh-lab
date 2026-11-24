/* holefill: close a boundary hole; one-shot fan for short convex loops,
 * concavity-driven block triangulation otherwise. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAXV 256
#define MAXF 512
#define MAXLOOP 64
#define MAXBLK 34

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { V3 v[MAXV]; int nv; Tri f[MAXF]; int nf; } Mesh;
typedef struct { int v[MAXLOOP]; int n; } Loop;
typedef struct { double x, y; } P2;

#define EPS 1e-9

typedef struct {
    int loop_n;
    int blocks;        /* sub-polygons used */
    int splits;        /* split vertices on the ring */
    int fill_tris;
    int self_int;      /* overlapping triangle pairs in the patch */
    int reflex;
    int k_from_concavity;
    int root;
    double concavity;
    double perimeter;
    int ok;
} FillStats;

/* ---------- mesh ---------- */

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

static int mesh_add_vertex(Mesh* m, float x, float y, float z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv].x = x; m->v[m->nv].y = y; m->v[m->nv].z = z;
    return m->nv++;
}

static int mesh_add_tri(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a < 0 || b < 0 || c < 0 || a >= m->nv || b >= m->nv || c >= m->nv) return -1;
    if (a == b || b == c || a == c) return -1;
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* ---------- 2D predicates ---------- */

static double orient2(P2 a, P2 b, P2 c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static int proper_cross(P2 a, P2 b, P2 c, P2 d) {
    double d1 = orient2(c, d, a), d2 = orient2(c, d, b);
    double d3 = orient2(a, b, c), d4 = orient2(a, b, d);
    return ((d1 > EPS && d2 < -EPS) || (d1 < -EPS && d2 > EPS)) &&
           ((d3 > EPS && d4 < -EPS) || (d3 < -EPS && d4 > EPS));
}

static int on_segment(P2 a, P2 b, P2 p) {
    if (fabs(orient2(a, b, p)) > EPS) return 0;
    double lo_x = a.x < b.x ? a.x : b.x, hi_x = a.x > b.x ? a.x : b.x;
    double lo_y = a.y < b.y ? a.y : b.y, hi_y = a.y > b.y ? a.y : b.y;
    return p.x >= lo_x - EPS && p.x <= hi_x + EPS &&
           p.y >= lo_y - EPS && p.y <= hi_y + EPS;
}

static int point_in_poly(P2 p, const P2* q, int n) {
    int inside = 0;
    for (int i = 0; i < n; i++) {
        P2 a = q[i], b = q[(i + 1) % n];
        if ((a.y > p.y) != (b.y > p.y)) {
            double x = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);
            if (x > p.x) inside = !inside;
        }
    }
    return inside;
}

static int point_strict_in_tri(P2 p, P2 a, P2 b, P2 c) {
    double d1 = orient2(a, b, p), d2 = orient2(b, c, p), d3 = orient2(c, a, p);
    return (d1 > EPS && d2 > EPS && d3 > EPS) ||
           (d1 < -EPS && d2 < -EPS && d3 < -EPS);
}

/* segment i-j usable as a split chord: no crossing, no touching a
 * non-incident vertex, midpoint inside the ring */
static int chord_ok(const P2* q, int n, int i, int j) {
    if (i == j) return 0;
    if ((i + 1) % n == j || (j + 1) % n == i) return 0;
    P2 a = q[i], b = q[j];
    double dx = b.x - a.x, dy = b.y - a.y;
    if (dx * dx + dy * dy < EPS * EPS) return 0;
    for (int k = 0; k < n; k++) {
        int k2 = (k + 1) % n;
        if (k == i || k == j || k2 == i || k2 == j) continue;
        if (proper_cross(a, b, q[k], q[k2])) return 0;
        if (on_segment(a, b, q[k]) || on_segment(a, b, q[k2])) return 0;
    }
    P2 mid = { (a.x + b.x) * 0.5, (a.y + b.y) * 0.5 };
    return point_in_poly(mid, q, n);
}

/* ---------- polygon metrics ---------- */

static double poly_area2(const P2* q, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        const P2* a = &q[i];
        const P2* b = &q[(i + 1) % n];
        s += a->x * b->y - b->x * a->y;
    }
    return s * 0.5;
}

static double poly_perimeter(const P2* q, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        const P2* a = &q[i];
        const P2* b = &q[(i + 1) % n];
        s += hypot(b->x - a->x, b->y - a->y);
    }
    return s;
}

static int cmp_pt(const void* pa, const void* pb) {
    const P2* a = (const P2*)pa;
    const P2* b = (const P2*)pb;
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    return 0;
}

static double hull_area(const P2* q, int n) {
    P2 p[MAXLOOP];
    memcpy(p, q, sizeof(P2) * (size_t)n);
    qsort(p, (size_t)n, sizeof(P2), cmp_pt);
    int k = 0;
    P2 h[2 * MAXLOOP];
    for (int i = 0; i < n; i++) {
        while (k >= 2 && orient2(h[k - 2], h[k - 1], p[i]) <= EPS) k--;
        h[k++] = p[i];
    }
    int lo = k;
    for (int i = n - 2; i >= 0; i--) {
        while (k > lo && orient2(h[k - 2], h[k - 1], p[i]) <= EPS) k--;
        h[k++] = p[i];
    }
    if (k < 4) return 0.0;
    return fabs(poly_area2(h, k - 1));
}

/* convexity deficit: 0 for a convex loop, ->1 as the loop collapses */
static double loop_concavity(const P2* q, int n) {
    double a = fabs(poly_area2(q, n));
    double h = hull_area(q, n);
    if (h <= EPS) return -1.0;
    double c = 1.0 - a / h;
    if (c < 0.0) c = 0.0;
    if (c > 1.0) c = 1.0;
    return c;
}

static int reflex_vertices(const P2* q, int n, int* out) {
    int r = 0;
    for (int i = 0; i < n; i++) {
        if (orient2(q[(i + n - 1) % n], q[i], q[(i + 1) % n]) < -EPS) out[r++] = i;
    }
    return r;
}

/* ---------- triangulation ---------- */

static int earclip(const P2* p, int m, int* out, int* nt) {
    int idx[2 * MAXLOOP];
    for (int i = 0; i < m; i++) idx[i] = i;
    if (poly_area2(p, m) < 0.0) {
        for (int i = 0; i < m; i++) idx[i] = m - 1 - i;
    }
    int cnt = m;
    int guard = 0;
    while (cnt > 3 && guard++ < 10000) {
        int found = 0;
        for (int t = 0; t < cnt; t++) {
            int i0 = idx[(t + cnt - 1) % cnt];
            int i1 = idx[t];
            int i2 = idx[(t + 1) % cnt];
            if (orient2(p[i0], p[i1], p[i2]) <= EPS) continue;
            int ok = 1;
            for (int k = 0; k < cnt; k++) {
                int kk = idx[k];
                if (kk == i0 || kk == i1 || kk == i2) continue;
                if (point_strict_in_tri(p[kk], p[i0], p[i1], p[i2])) { ok = 0; break; }
            }
            if (!ok) continue;
            out[3 * (*nt) + 0] = i0; out[3 * (*nt) + 1] = i1; out[3 * (*nt) + 2] = i2;
            (*nt)++;
            for (int k = t; k < cnt - 1; k++) idx[k] = idx[k + 1];
            cnt--;
            found = 1;
            break;
        }
        if (!found) return -1;
    }
    if (cnt != 3) return -1;
    if (orient2(p[idx[0]], p[idx[1]], p[idx[2]]) <= EPS) return -1;
    out[3 * (*nt) + 0] = idx[0]; out[3 * (*nt) + 1] = idx[1]; out[3 * (*nt) + 2] = idx[2];
    (*nt)++;
    return 0;
}

static int tri_overlap(P2 a0, P2 a1, P2 a2, P2 b0, P2 b1, P2 b2) {
    P2 A[3] = { a0, a1, a2 };
    P2 B[3] = { b0, b1, b2 };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (proper_cross(A[i], A[(i + 1) % 3], B[j], B[(j + 1) % 3])) return 1;
        }
    }
    for (int i = 0; i < 3; i++) {
        if (point_strict_in_tri(A[i], b0, b1, b2)) return 1;
        if (point_strict_in_tri(B[i], a0, a1, a2)) return 1;
    }
    return 0;
}

/* ---------- boundary loops ---------- */

typedef struct { int a, b, count; } Edge;

static int edge_find(Edge* e, int ne, int a, int b) {
    if (a > b) { int t = a; a = b; b = t; }
    for (int i = 0; i < ne; i++) if (e[i].a == a && e[i].b == b) return i;
    return -1;
}

static int mesh_boundary_loops(const Mesh* m, Loop* loops, int max_loops) {
    Edge e[3 * MAXF];
    int ne = 0;
    for (int i = 0; i < m->nf; i++) {
        int v[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
        for (int k = 0; k < 3; k++) {
            int a = v[k], b = v[(k + 1) % 3];
            int id = edge_find(e, ne, a, b);
            if (id < 0) {
                if (ne >= 3 * MAXF) return -1;
                int lo = a < b ? a : b, hi = a < b ? b : a;
                e[ne].a = lo; e[ne].b = hi; e[ne].count = 1; ne++;
            } else {
                e[id].count++;
            }
        }
    }
    int ha[3 * MAXF], hb[3 * MAXF], nh = 0;
    for (int i = 0; i < m->nf; i++) {
        int v[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
        for (int k = 0; k < 3; k++) {
            int a = v[k], b = v[(k + 1) % 3];
            int id = edge_find(e, ne, a, b);
            if (id >= 0 && e[id].count == 1) { ha[nh] = a; hb[nh] = b; nh++; }
        }
    }
    int used[3 * MAXF];
    memset(used, 0, sizeof(used));
    int nl = 0;
    for (int s = 0; s < nh; s++) {
        if (used[s] || nl >= max_loops) continue;
        Loop* L = &loops[nl];
        L->n = 0;
        int cur = s;
        while (1) {
            used[cur] = 1;
            if (L->n >= MAXLOOP) return -1;
            L->v[L->n++] = ha[cur];
            int next = -1;
            for (int k = 0; k < nh; k++) {
                if (ha[k] == hb[cur]) { next = k; break; }
            }
            if (next < 0) return -1;      /* dangling boundary */
            if (next == s) break;         /* closed the loop */
            if (used[next]) return -1;    /* not a simple cycle */
            cur = next;
        }
        if (L->n < 3) return -1;
        nl++;
    }
    return nl;
}

/* ---------- fill ---------- */

/* count how many vertices the root can see across the interior */
static int visibility_count(const P2* q, int n, int r) {
    int c = 0;
    for (int j = 0; j < n; j++) if (chord_ok(q, n, r, j)) c++;
    return c;
}

/* candidate roots, most visible first; reflex vertices preferred */
static int rank_roots(const P2* q, int n, int* roots) {
    int refl[MAXLOOP];
    int nr = reflex_vertices(q, n, refl);
    int nc = 0;
    if (nr == 0) {
        for (int i = 0; i < n; i++) roots[nc++] = i;
    } else {
        for (int i = 0; i < nr; i++) roots[nc++] = refl[i];
    }
    int vis[MAXLOOP];
    for (int i = 0; i < nc; i++) vis[i] = visibility_count(q, n, roots[i]);
    for (int i = 0; i < nc; i++) {
        int best = i;
        for (int j = i + 1; j < nc; j++)
            if (vis[j] > vis[best] || (vis[j] == vis[best] && roots[j] < roots[best])) best = j;
        int t = roots[i]; roots[i] = roots[best]; roots[best] = t;
        t = vis[i]; vis[i] = vis[best]; vis[best] = t;
    }
    return nc;
}

/* roots that give no split at all are useless; raise k by one until it fits */
static int place_splits(const P2* q, int n, int root, int k, int* pos, int* sel) {
    if (k <= 1) { pos[0] = root; return 0; }
    double d[MAXLOOP + 1];
    int idx[MAXLOOP + 1];
    idx[0] = root;
    d[0] = 0.0;
    for (int i = 0; i < n; i++) {
        int a = (root + i) % n, b = (root + i + 1) % n;
        idx[i + 1] = b;
        d[i + 1] = d[i] + hypot(q[b].x - q[a].x, q[b].y - q[a].y);
    }
    double P = d[n];
    int got = 0;
    for (int j = 1; j < k; j++) {
        double target = P * (double)j / (double)k;
        int lo = got ? sel[got - 1] + 1 : 1;
        int best = -1;
        double bd = 1e30;
        for (int c = lo; c < n; c++) {
            if (n - c < k - 1 - j) continue;
            double dv = fabs(d[c] - target);
            if (dv >= bd) continue;
            if (!chord_ok(q, n, root, idx[c])) continue;
            bd = dv; best = c;
        }
        if (best < 0) return -1;
        sel[got++] = best;
    }
    for (int i = 0; i < n + 1; i++) pos[i] = idx[i];
    return got;
}

/* split the ring into k pieces around root; each piece is root + an arc run.
 * returns piece count, or -1 if the placement does not work */
static int make_pieces(int n, const int* idx, const int* sel, int m,
                       int pieces[][MAXLOOP], int* plen) {
    int bounds[MAXLOOP];
    bounds[0] = 0;
    for (int i = 0; i < m; i++) bounds[i + 1] = sel[i];
    bounds[m + 1] = n;
    for (int t = 0; t <= m; t++) {
        int len = 0;
        if (t > 0 && t < m) { pieces[t][0] = idx[0]; len = 1; }
        for (int p = bounds[t]; p <= bounds[t + 1]; p++) {
            int v = idx[p];
            if (len > 0 && v == pieces[t][len - 1]) continue;
            if (len >= MAXLOOP) return -1;
            pieces[t][len++] = v;
        }
        if (len < 3) return -1;
        plen[t] = len;
    }
    return m + 1;
}

static int fill_with_k(const P2* q, int n, int k, int* tris, int* nt,
                       int* used_pieces, int* used_root) {
    int roots[MAXLOOP];
    int nr = rank_roots(q, n, roots);
    for (int ri = 0; ri < nr; ri++) {
        int idx[MAXLOOP + 1];
        int sel[MAXLOOP];
        int m = place_splits(q, n, roots[ri], k, idx, sel);
        if (m < 0) continue;
        int pieces[MAXLOOP][MAXLOOP], plen[MAXLOOP];
        int nb = make_pieces(n, idx, sel, m, pieces, plen);
        if (nb < 0) continue;
        int local[3 * MAXLOOP], lt_all = 0;
        int ok = 1;
        int all[3 * MAXLOOP];
        for (int b = 0; b < nb; b++) {
            P2 sub[MAXLOOP];
            int verts[MAXLOOP];
            for (int i = 0; i < plen[b]; i++) verts[i] = pieces[b][i];
            for (int i = 0; i < plen[b]; i++) sub[i] = q[verts[i]];
            if (poly_area2(sub, plen[b]) < 0.0) {
                for (int i = 0; i < plen[b] / 2; i++) {
                    P2 t2 = sub[i]; sub[i] = sub[plen[b] - 1 - i]; sub[plen[b] - 1 - i] = t2;
                    int t3 = verts[i]; verts[i] = verts[plen[b] - 1 - i]; verts[plen[b] - 1 - i] = t3;
                }
            }
            int lt = 0;
            if (earclip(sub, plen[b], local, &lt) != 0) { ok = 0; break; }
            for (int i = 0; i < lt; i++) {
                all[3 * lt_all + 0] = verts[local[3 * i + 0]];
                all[3 * lt_all + 1] = verts[local[3 * i + 1]];
                all[3 * lt_all + 2] = verts[local[3 * i + 2]];
                lt_all++;
            }
        }
        if (!ok) continue;
        memcpy(tris, all, sizeof(int) * 3 * (size_t)lt_all);
        *nt = lt_all;
        *used_pieces = nb;
        *used_root = roots[ri];
        return 0;
    }
    return -1;
}

static int count_self_int(const P2* q, const int* tris, int nt) {
    int c = 0;
    for (int i = 0; i < nt; i++) {
        for (int j = i + 1; j < nt; j++) {
            if (tri_overlap(q[tris[3 * i + 0]], q[tris[3 * i + 1]], q[tris[3 * i + 2]],
                            q[tris[3 * j + 0]], q[tris[3 * j + 1]], q[tris[3 * j + 2]])) c++;
        }
    }
    return c;
}

static void loop_to_plane(const Mesh* mesh, const Loop* L, P2* q, int* pos) {
    /* Newell normal */
    double nx = 0.0, ny = 0.0, nz = 0.0;
    for (int i = 0; i < L->n; i++) {
        V3 p = mesh->v[L->v[i]];
        V3 r = mesh->v[L->v[(i + 1) % L->n]];
        nx += ((double)p.y - r.y) * ((double)p.z + r.z);
        ny += ((double)p.z - r.z) * ((double)p.x + r.x);
        nz += ((double)p.x - r.x) * ((double)p.y + r.y);
    }
    double len = sqrt(nx * nx + ny * ny + nz * nz);
    if (len < 1e-12) len = 1.0;
    nx /= len; ny /= len; nz /= len;
    double ax = fabs(nx), ay = fabs(ny), az = fabs(nz);
    double ux, uy, uz;
    if (ax <= ay && ax <= az) { ux = 0; uy = -nz; uz = ny; }
    else if (ay <= az)        { ux = nz; uy = 0; uz = -nx; }
    else                      { ux = -ny; uy = nx; uz = 0; }
    double ul = sqrt(ux * ux + uy * uy + uz * uz);
    if (ul < 1e-12) ul = 1.0;
    ux /= ul; uy /= ul; uz /= ul;
    double wx = ny * uz - nz * uy;
    double wy = nz * ux - nx * uz;
    double wz = nx * uy - ny * ux;
    for (int i = 0; i < L->n; i++) {
        V3 p = mesh->v[L->v[i]];
        q[i].x = (double)p.x * ux + (double)p.y * uy + (double)p.z * uz;
        q[i].y = (double)p.x * wx + (double)p.y * wy + (double)p.z * wz;
        pos[i] = i;
    }
}

/* fill one boundary loop into the mesh; stats -> st */
static int hole_fill(Mesh* mesh, const Loop* L, FillStats* st) {
    memset(st, 0, sizeof(*st));
    st->loop_n = L->n;
    if (L->n < 3 || L->n > MAXLOOP) return -1;

    P2 q[MAXLOOP];
    int pos[MAXLOOP];
    loop_to_plane(mesh, L, q, pos);

    /* make the 2D copy CCW; remember whether that reversed the loop */
    int reversed = 0;
    if (poly_area2(q, L->n) < 0.0) {
        for (int i = 0; i < L->n / 2; i++) {
            P2 t = q[i]; q[i] = q[L->n - 1 - i]; q[L->n - 1 - i] = t;
            int s = pos[i]; pos[i] = pos[L->n - 1 - i]; pos[L->n - 1 - i] = s;
        }
        reversed = 1;
    }

    st->concavity = loop_concavity(q, L->n);
    st->perimeter = poly_perimeter(q, L->n);
    if (st->concavity < 0.0) return -1;

    int refl[MAXLOOP];
    int nr = reflex_vertices(q, L->n, refl);
    st->reflex = nr;

    /* chunks scale with the concavity; at least two once the loop is concave,
     * and never so many that a chunk would drop below three vertices */
    int k = 1;
    if (nr > 0) {
        k = 1 + (int)floor(st->concavity / 0.08 + 0.5);
        if (k < 2) k = 2;
        int hi = L->n / 2;
        if (k > hi) k = hi;
        if (k < 2) k = 2;
    }

    int tris[3 * MAXLOOP];      /* any patch of an n-gon has exactly n-2 triangles */
    int nt = 0, nb = 0, root = -1;
    int used_k = 1;
    int ok = 0;
    if (k <= 1) {
        /* short and convex: single one-shot fan */
        for (int i = 1; i < L->n - 1; i++) {
            tris[3 * nt + 0] = 0; tris[3 * nt + 1] = i; tris[3 * nt + 2] = i + 1; nt++;
        }
        nb = 1; ok = 1;
    } else {
        for (int kk = k; kk >= 2 && !ok; kk--) {
            if (fill_with_k(q, L->n, kk, tris, &nt, &nb, &root) == 0) { used_k = kk; ok = 1; }
        }
        if (!ok) {
            /* degenerate layout: fall back to one ear-clipped patch */
            P2 sub[MAXLOOP];
            int verts[MAXLOOP];
            for (int i = 0; i < L->n; i++) { sub[i] = q[i]; verts[i] = i; }
            int local[3 * MAXLOOP], lt = 0;
            if (earclip(sub, L->n, local, &lt) == 0) {
                for (int i = 0; i < lt; i++) {
                    tris[3 * nt + 0] = verts[local[3 * i + 0]];
                    tris[3 * nt + 1] = verts[local[3 * i + 1]];
                    tris[3 * nt + 2] = verts[local[3 * i + 2]];
                    nt++;
                }
                nb = 1; used_k = 1; ok = 1;
            }
        }
    }
    if (!ok) return -1;

    st->blocks = nb;
    st->splits = nb > 1 ? nb - 1 : 0;
    st->fill_tris = nt;
    st->self_int = count_self_int(q, tris, nt);
    st->k_from_concavity = used_k;
    st->root = root;
    st->ok = 1;

    /* append with the winding that opposes the boundary half-edges */
    for (int i = 0; i < nt; i++) {
        int a = L->v[pos[tris[3 * i + 0]]];
        int b = L->v[pos[tris[3 * i + 1]]];
        int c = L->v[pos[tris[3 * i + 2]]];
        int r;
        if (reversed) r = mesh_add_tri(mesh, a, b, c);
        else          r = mesh_add_tri(mesh, a, c, b);
        if (r < 0) return -1;
    }
    return 0;
}

/* ---------- edge statistics ---------- */

static void mesh_edge_stats(const Mesh* mesh, int* open, int* nonman, int* non_two) {
    Edge e[3 * MAXF];
    int ne = 0;
    for (int i = 0; i < mesh->nf; i++) {
        int v[3] = { mesh->f[i].a, mesh->f[i].b, mesh->f[i].c };
        for (int k = 0; k < 3; k++) {
            int a = v[k], b = v[(k + 1) % 3];
            int id = edge_find(e, ne, a, b);
            if (id < 0) {
                int lo = a < b ? a : b, hi = a < b ? b : a;
                e[ne].a = lo; e[ne].b = hi; e[ne].count = 1; ne++;
            } else {
                e[id].count++;
            }
        }
    }
    int o = 0, nm = 0, nt2 = 0;
    for (int i = 0; i < ne; i++) {
        if (e[i].count == 1) o++;
        if (e[i].count > 2) nm++;
        if (e[i].count != 2) nt2++;
    }
    *open = o; *nonman = nm; *non_two = nt2;
}

/* each edge must be used once in each direction */
static int mesh_oriented(const Mesh* mesh) {
    Edge e[6 * MAXF];
    int ne = 0;
    for (int i = 0; i < mesh->nf; i++) {
        int v[3] = { mesh->f[i].a, mesh->f[i].b, mesh->f[i].c };
        for (int k = 0; k < 3; k++) {
            int a = v[k], b = v[(k + 1) % 3];
            int lo = a < b ? a : b, hi = a < b ? b : a;
            int dir = (a < b) ? 1 : -1;
            int id = -1;
            for (int t = 0; t < ne; t++) {
                if (e[t].a == lo && e[t].b == hi) { id = t; break; }
            }
            if (id < 0) { e[ne].a = lo; e[ne].b = hi; e[ne].count = dir; ne++; }
            else e[id].count += dir;
        }
    }
    for (int i = 0; i < ne; i++) if (e[i].count != 0) return 0;
    return 1;
}

/* ---------- fixtures ---------- */

#define OUTER_SCALE 3.0

/* open-bottomed shell: hole ring on z=0, outer ring, wall, top fan */
static int build_shell(Mesh* m, const P2* ring, int n, float zjitter) {
    mesh_init(m);
    for (int i = 0; i < n; i++) {
        float z = 0.0f;
        if (zjitter != 0.0f) z = (i & 1) ? zjitter : -zjitter;
        mesh_add_vertex(m, (float)ring[i].x, (float)ring[i].y, z);
    }
    for (int i = 0; i < n; i++)
        mesh_add_vertex(m, (float)(ring[i].x * OUTER_SCALE), (float)(ring[i].y * OUTER_SCALE), 0.0f);
    for (int i = 0; i < n; i++)
        mesh_add_vertex(m, (float)(ring[i].x * OUTER_SCALE), (float)(ring[i].y * OUTER_SCALE), 1.0f);
    int apex = mesh_add_vertex(m, 0.0f, 0.0f, 1.8f);
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        mesh_add_tri(m, n + i, i, j);              /* bottom annulus */
        mesh_add_tri(m, n + i, j, n + j);
        mesh_add_tri(m, n + i, n + j, 2 * n + j);  /* wall */
        mesh_add_tri(m, n + i, 2 * n + j, 2 * n + i);
        mesh_add_tri(m, 2 * n + i, 2 * n + j, apex); /* top fan */
    }
    return 0;
}

/* unit cube with top and bottom faces removed -> two 4-vertex holes */
static int build_open_box(Mesh* m) {
    mesh_init(m);
    mesh_add_vertex(m, 0, 0, 0); mesh_add_vertex(m, 1, 0, 0);
    mesh_add_vertex(m, 1, 1, 0); mesh_add_vertex(m, 0, 1, 0);
    mesh_add_vertex(m, 0, 0, 1); mesh_add_vertex(m, 1, 0, 1);
    mesh_add_vertex(m, 1, 1, 1); mesh_add_vertex(m, 0, 1, 1);
    int f[][3] = {
        {0,1,5},{0,5,4},      /* -y */
        {1,2,6},{1,6,5},      /* +x */
        {2,3,7},{2,7,6},      /* +y */
        {3,0,4},{3,4,7},      /* -x */
        /* +z and -z omitted */
    };
    for (int i = 0; i < (int)(sizeof(f)/sizeof(f[0])); i++)
        mesh_add_tri(m, f[i][0], f[i][1], f[i][2]);
    return 0;
}

/* ---------- tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void ring_regular(P2* q, int n, double r) {
    for (int i = 0; i < n; i++) {
        double a = 2.0 * M_PI * i / n;
        q[i].x = r * cos(a); q[i].y = r * sin(a);
    }
}

static void ring_medium(P2* q) {
    double r[12] = { 1.20, 1.20, 0.62, 1.20, 1.20, 0.62, 1.20, 1.20, 0.62, 1.20, 1.20, 0.62 };
    for (int i = 0; i < 12; i++) {
        double a = 2.0 * M_PI * i / 12.0;
        q[i].x = r[i] * cos(a); q[i].y = r[i] * sin(a);
    }
}

static void ring_large(P2* q) {
    double r[16];
    for (int i = 0; i < 16; i++) r[i] = (i % 2 == 0) ? 2.00 : 0.85;
    for (int i = 0; i < 16; i++) {
        double a = 2.0 * M_PI * i / 16.0;
        q[i].x = r[i] * cos(a); q[i].y = r[i] * sin(a);
    }
}

/* the naive rule under test: chunk count from perimeter, splits by arc length */
static void control_perimeter(Mesh* m, const Loop* L, int* tris, double* self_int,
                              double* conc, double* perim, int* kp, int* bad_blocks) {
    P2 q[MAXLOOP];
    int pos[MAXLOOP];
    loop_to_plane(m, L, q, pos);
    if (poly_area2(q, L->n) < 0.0) {
        for (int i = 0; i < L->n / 2; i++) { P2 t = q[i]; q[i] = q[L->n-1-i]; q[L->n-1-i] = t; }
    }
    int n = L->n;
    *conc = loop_concavity(q, n);
    *perim = poly_perimeter(q, n);
    int k = (int)ceil(*perim / 4.0);
    if (k < 1) k = 1;
    *kp = k;
    double d[MAXLOOP + 1];
    d[0] = 0.0;
    for (int i = 0; i < n; i++) d[i + 1] = d[i] + hypot(q[(i + 1) % n].x - q[i].x, q[(i + 1) % n].y - q[i].y);
    int sel[MAXLOOP], m_sel = 0;
    if (k >= 2) {
        for (int j = 1; j < k; j++) {
            double target = *perim * (double)j / (double)k;
            int best = -1; double bd = 1e30;
            for (int v = 0; v < n; v++) {
                double dv = fabs(d[v] - target);
                if (dv >= bd) continue;
                int clash = 0;
                for (int s = 0; s < m_sel; s++) {
                    int dd = v - sel[s]; if (dd < 0) dd = -dd;
                    if (dd > n - dd) dd = n - dd;
                    if (dd < 2) { clash = 1; break; }
                }
                if (clash) continue;
                bd = dv; best = v;
            }
            if (best >= 0) sel[m_sel++] = best;
        }
    }
    for (int i = 0; i < m_sel; i++)
        for (int j = i + 1; j < m_sel; j++)
            if (sel[j] < sel[i]) { int t = sel[i]; sel[i] = sel[j]; sel[j] = t; }
    int blocks[MAXBLK][MAXLOOP], blen[MAXBLK], nb = 0;
    if (m_sel == 0) {
        for (int i = 0; i < n; i++) blocks[0][i] = i;
        blen[0] = n; nb = 1;
    } else {
        for (int t = 0; t <= m_sel; t++) {
            int lo = (t == 0) ? 0 : sel[t - 1];
            int hi = (t == m_sel) ? n : sel[t];
            int len = 0;
            for (int p = lo; p <= hi; p++) blocks[nb][len++] = p % n;
            if (len < 3) continue;
            blen[nb] = len; nb++;
        }
        if (nb == 0) { *self_int = -1; *tris = 0; *bad_blocks = 0; return; }
    }
    int out[3 * MAXLOOP], nt = 0, bad = 0;
    for (int b = 0; b < nb; b++) {
        /* one-shot fan per chunk: no interior check, crossing blocks stay crossing */
        if (blen[b] >= 3) {
            for (int i = 1; i < blen[b] - 1; i++) {
                out[3 * nt + 0] = blocks[b][0];
                out[3 * nt + 1] = blocks[b][i];
                out[3 * nt + 2] = blocks[b][i + 1];
                nt++;
            }
        } else {
            bad++;
        }
    }
    *tris = nt;
    *bad_blocks = bad;
    *self_int = count_self_int(q, out, nt);
}

int main(void) {
    P2 small[6], medium[12], large[16];
    ring_regular(small, 6, 0.55);
    ring_medium(medium);
    ring_large(large);

    /* ---- geometry/law checks on the projected loops ---- */
    check("small_concavity_zero", fabs(loop_concavity(small, 6)) < 1e-6);
    check("medium_concavity", fabs(loop_concavity(medium, 12) - 0.255748) < 1e-4);
    check("large_concavity", fabs(loop_concavity(large, 16) - 0.539983) < 1e-4);

    struct { const char* name; P2* q; int n; float zj; int exp_blocks; int exp_tris;
             double exp_si; double exp_c; } cases[] = {
        { "small",  small,   6, 0.0f, 1,  4, 0.0, 0.000000 },
        { "medium", medium, 12, 0.0f, 4, 10, 0.0, 0.255748 },
        { "large",  large,  16, 0.0f, 8, 14, 0.0, 0.539983 },
    };

    for (int ci = 0; ci < 3; ci++) {
        Mesh m;
        build_shell(&m, cases[ci].q, cases[ci].n, cases[ci].zj);
        Loop loops[4];
        int nl = mesh_boundary_loops(&m, loops, 4);
        int open0, nm0, nt0;
        mesh_edge_stats(&m, &open0, &nm0, &nt0);
        char tag[64];

        snprintf(tag, sizeof(tag), "%s_one_loop", cases[ci].name);
        check(tag, nl == 1 && loops[0].n == cases[ci].n);
        snprintf(tag, sizeof(tag), "%s_open_before", cases[ci].name);
        check(tag, open0 == cases[ci].n && nm0 == 0);

        FillStats st;
        int rc = hole_fill(&m, &loops[0], &st);
        snprintf(tag, sizeof(tag), "%s_fill_ok", cases[ci].name);
        check(tag, rc == 0 && st.ok);

        int open1, nm1, nt1;
        mesh_edge_stats(&m, &open1, &nm1, &nt1);
        snprintf(tag, sizeof(tag), "%s_closed", cases[ci].name);
        check(tag, open1 == 0);
        snprintf(tag, sizeof(tag), "%s_no_nonmanifold", cases[ci].name);
        check(tag, nm1 == 0 && nt1 == 0);
        snprintf(tag, sizeof(tag), "%s_oriented", cases[ci].name);
        check(tag, mesh_oriented(&m));
        snprintf(tag, sizeof(tag), "%s_tri_count", cases[ci].name);
        check(tag, st.fill_tris == cases[ci].exp_tris);
        snprintf(tag, sizeof(tag), "%s_blocks", cases[ci].name);
        check(tag, st.blocks == cases[ci].exp_blocks);
        snprintf(tag, sizeof(tag), "%s_self_int", cases[ci].name);
        check(tag, st.self_int == (int)cases[ci].exp_si);
        snprintf(tag, sizeof(tag), "%s_concavity", cases[ci].name);
        check(tag, fabs(st.concavity - cases[ci].exp_c) < 1e-4);
        snprintf(tag, sizeof(tag), "%s_chunk_count", cases[ci].name);
        check(tag, st.k_from_concavity == cases[ci].exp_blocks);
        snprintf(tag, sizeof(tag), "%s_no_new_vertices", cases[ci].name);
        check(tag, m.nv == 3 * cases[ci].n + 1);

        printf("%s_mesh: fill_tris=%d self_int=%d nonmanifold=%d blocks=%d splits=%d concavity=%.6f reflex=%d\n",
               cases[ci].name, st.fill_tris, st.self_int, nm1, st.blocks, st.splits,
               st.concavity, st.reflex);
    }

    /* ---- small loop with one shallow notch: still one piece, must not cross ---- */
    {
        P2 q[7] = {
            { 0.6172, 0.0000 }, { 0.2917, 0.3658 }, { -0.1044, 0.4575 },
            { -0.7057, 0.3399 }, { -0.7674, -0.3696 }, { -0.2339, -1.0247 },
            { 0.2660, -0.3336 }
        };
        Mesh m;
        build_shell(&m, q, 7, 0.0f);
        Loop loops[4];
        int nl = mesh_boundary_loops(&m, loops, 4);
        FillStats st;
        int rc = hole_fill(&m, &loops[0], &st);
        int open1, nm1, nt1;
        mesh_edge_stats(&m, &open1, &nm1, &nt1);
        check("shallow_notch_one_loop", nl == 1 && loops[0].n == 7);
        check("shallow_notch_fill_ok", rc == 0 && st.ok);
        check("shallow_notch_closed", open1 == 0 && nt1 == 0);
        check("shallow_notch_self_int", st.self_int == 0);
        check("shallow_notch_tris", st.fill_tris == 5);
        printf("shallow_notch: concavity=%.6f reflex=%d chunks=%d root=%d self_int=%d\n",
               st.concavity, st.reflex, st.k_from_concavity, st.root, st.self_int);
    }

    /* ---- non-planar boundary: same loop with alternating z jitter ---- */
    {
        Mesh m;
        build_shell(&m, large, 16, 0.03f);
        Loop loops[4];
        int nl = mesh_boundary_loops(&m, loops, 4);
        check("jitter_one_loop", nl == 1 && loops[0].n == 16);
        FillStats st;
        int rc = hole_fill(&m, &loops[0], &st);
        int open1, nm1, nt1;
        mesh_edge_stats(&m, &open1, &nm1, &nt1);
        check("jitter_fill_ok", rc == 0 && st.ok);
        check("jitter_closed", open1 == 0 && nt1 == 0);
        check("jitter_self_int", st.self_int == 0);
        check("jitter_tris", st.fill_tris == 14);
        printf("jitter_mesh: fill_tris=%d self_int=%d nonmanifold=%d\n", st.fill_tris, st.self_int, nm1);
    }

    /* ---- two holes in one mesh: find both, fill both, close ---- */
    {
        Mesh m;
        build_open_box(&m);
        Loop loops[4];
        int nl = mesh_boundary_loops(&m, loops, 4);
        check("box_two_loops", nl == 2 && loops[0].n == 4 && loops[1].n == 4);
        FillStats st;
        int rc1 = hole_fill(&m, &loops[0], &st);
        check("box_fill_first", rc1 == 0 && st.fill_tris == 2);
        nl = mesh_boundary_loops(&m, loops, 4);
        check("box_one_loop_left", nl == 1 && loops[0].n == 4);
        int rc2 = hole_fill(&m, &loops[0], &st);
        check("box_fill_second", rc2 == 0 && st.fill_tris == 2);
        int open1, nm1, nt1;
        mesh_edge_stats(&m, &open1, &nm1, &nt1);
        check("box_closed", open1 == 0 && nm1 == 0 && nt1 == 0);
        check("box_oriented", mesh_oriented(&m));
    }

    /* ---- counterfactual: chunk count from perimeter proportion ---- */
    {
        double self[3], conc[3], per[3];
        int tris[3], kp[3], bad[3];
        const char* names[3] = { "small", "medium", "large" };
        P2* rings[3] = { small, medium, large };
        int ns[3] = { 6, 12, 16 };
        for (int i = 0; i < 3; i++) {
            Mesh m;
            build_shell(&m, rings[i], ns[i], 0.0f);
            Loop loops[4];
            mesh_boundary_loops(&m, loops, 4);
            control_perimeter(&m, &loops[0], &tris[i], &self[i], &conc[i], &per[i], &kp[i], &bad[i]);
        }
        check("control_small_self_int", self[0] == 0);
        check("control_medium_self_int_gt0", self[1] > 0);
        check("control_large_self_int_gt0", self[2] > 0);
        for (int i = 0; i < 3; i++)
            printf("control_%s: P=%.6f concavity=%.6f k_perimeter=%d fill_tris=%d self_int=%.0f bad_blocks=%d\n",
                   names[i], per[i], conc[i], kp[i], tris[i], self[i], bad[i]);
    }

    /* ---- summary: the nine numbers, with the chunk counts side by side ---- */
    {
        const char* names[3] = { "small", "medium", "large" };
        P2* rings[3] = { small, medium, large };
        int ns[3] = { 6, 12, 16 };
        for (int i = 0; i < 3; i++) {
            Mesh m;
            build_shell(&m, rings[i], ns[i], 0.0f);
            Loop loops[4];
            mesh_boundary_loops(&m, loops, 4);
            FillStats st;
            hole_fill(&m, &loops[0], &st);
            int open1, nm1, nt1;
            mesh_edge_stats(&m, &open1, &nm1, &nt1);
            int kp = (int)ceil(st.perimeter / 4.0);
            if (kp < 1) kp = 1;
            printf("summary_%s: fill_tris=%d self_int=%d nonmanifold_edges=%d "
                   "concavity=%.6f k_concavity=%d k_perimeter=%d reflex=%d\n",
                   names[i], st.fill_tris, st.self_int, nt1,
                   st.concavity, st.k_from_concavity, kp, st.reflex);
        }
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
