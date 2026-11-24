
/* normal: area-weighted vertex normals, full and incremental grouped updates. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 32
#define MAXF 32

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

static V3 sub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 add3(V3 p, V3 q) { V3 r = { p.x+q.x, p.y+q.y, p.z+q.z }; return r; }
static V3 scale3(V3 p, float s) { V3 r = { p.x*s, p.y*s, p.z*s }; return r; }
static V3 cross(V3 p, V3 q) { V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r; }
static float dot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float len(V3 p) { return sqrtf(dot(p,p)); }
static V3 norm(V3 p) { float l = len(p); V3 r = {0,0,0}; if (l > 0) { r.x=p.x/l; r.y=p.y/l; r.z=p.z/l; } return r; }

typedef struct { V3 v[MAXV]; int nv; Tri f[MAXF]; int nf; } Mesh;

static int add_v(Mesh* m, float x, float y, float z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv].x = x; m->v[m->nv].y = y; m->v[m->nv].z = z; return m->nv++;
}
static int add_f(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a<0||b<0||c<0||a>=m->nv||b>=m->nv||c>=m->nv) return -1;
    m->f[m->nf].a=a; m->f[m->nf].b=b; m->f[m->nf].c=c; return m->nf++;
}

/* area-weighted vertex normals: accumulate un-normalized face normals */
static int vertex_normals(const Mesh* m, V3* out) {
    if (m->nf == 0) return -1;
    for (int i = 0; i < m->nv; i++) { out[i].x = out[i].y = out[i].z = 0.0f; }
    for (int i = 0; i < m->nf; i++) {
        V3 a = m->v[m->f[i].a], b = m->v[m->f[i].b], c = m->v[m->f[i].c];
        V3 n = cross(sub(b, a), sub(c, a));   /* length = 2*area */
        out[m->f[i].a].x += n.x; out[m->f[i].a].y += n.y; out[m->f[i].a].z += n.z;
        out[m->f[i].b].x += n.x; out[m->f[i].b].y += n.y; out[m->f[i].b].z += n.z;
        out[m->f[i].c].x += n.x; out[m->f[i].c].y += n.y; out[m->f[i].c].z += n.z;
    }
    int bad = 0;
    for (int i = 0; i < m->nv; i++) {
        float l = len(out[i]);
        if (l < 1e-12f) { bad++; continue; }
        out[i] = norm(out[i]);
    }
    return bad;   /* number of vertices with zero-length accumulator */
}

/*
 * Incremental cache.  accum[i] is the unnormalized sum of face normals at
 * vertex index i.  Several indices may share one normal group; their accum
 * vectors are summed before normalization and every group member receives the
 * same normal.  This models welded/co-located normal corners.
 *
 * group_eps < 0 gives singleton groups and therefore the old per-index result.
 * group_eps >= 0 puts co-located vertices (distance <= eps) in one group.
 */
typedef struct {
    Mesh mesh;
    int incident[MAXV][MAXF];
    int incident_n[MAXV];
    V3 face_n[MAXF];
    V3 accum[MAXV];
    V3 normal[MAXV];

    int group[MAXV];
    unsigned members[MAXV];
    int position_groups;
    float group_eps;

    unsigned dirty_mask;
    int updated_faces;
} NormalCache;

static int dsu_find(int parent[MAXV], int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static void dsu_unite(int parent[MAXV], int a, int b) {
    a = dsu_find(parent, a);
    b = dsu_find(parent, b);
    if (a != b) parent[b] = a;
}

static void build_position_groups(const Mesh* m, float eps,
                                  int group[MAXV], unsigned members[MAXV]) {
    int parent[MAXV];
    for (int i = 0; i < MAXV; i++) parent[i] = i;

    if (eps >= 0.0f) {
        float eps2 = eps * eps;
        for (int i = 0; i < m->nv; i++) {
            for (int j = i + 1; j < m->nv; j++) {
                float dx = m->v[i].x - m->v[j].x;
                float dy = m->v[i].y - m->v[j].y;
                float dz = m->v[i].z - m->v[j].z;
                if (dx*dx + dy*dy + dz*dz <= eps2) dsu_unite(parent, i, j);
            }
        }
    }

    for (int i = 0; i < m->nv; i++) group[i] = dsu_find(parent, i);
    for (int i = 0; i < MAXV; i++) members[i] = 0;
    for (int i = 0; i < m->nv; i++) members[group[i]] |= 1u << i;
}

/* Full reference for grouped normals, used by the differential self-tests. */
static int grouped_vertex_normals(const Mesh* m, const int* group, V3* out) {
    V3 accum[MAXV] = {0};
    V3 group_sum[MAXV] = {0};

    for (int f = 0; f < m->nf; f++) {
        Tri t = m->f[f];
        V3 n = cross(sub(m->v[t.b], m->v[t.a]), sub(m->v[t.c], m->v[t.a]));
        accum[t.a] = add3(accum[t.a], n);
        accum[t.b] = add3(accum[t.b], n);
        accum[t.c] = add3(accum[t.c], n);
    }
    for (int i = 0; i < m->nv; i++) {
        group_sum[group[i]] = add3(group_sum[group[i]], accum[i]);
    }

    int bad = 0;
    for (int i = 0; i < m->nv; i++) {
        float l = len(group_sum[group[i]]);
        if (l < 1e-12f) {
            out[i] = (V3){0,0,0};
            bad++;
        } else {
            out[i] = scale3(group_sum[group[i]], 1.0f / l);
        }
    }
    return bad;
}

static void cache_recompute_group(NormalCache* c, int g) {
    V3 sum = {0,0,0};
    unsigned mask = c->members[g];
    for (int i = 0; i < c->mesh.nv; i++) {
        if (mask & (1u << i)) sum = add3(sum, c->accum[i]);
    }

    float l = len(sum);
    V3 n = {0,0,0};
    if (!(l < 1e-12f)) n = norm(sum);  /* norm() also guards l == 0 */

    for (int i = 0; i < c->mesh.nv; i++) {
        if (mask & (1u << i)) c->normal[i] = n;
    }
}

static int normal_cache_init(NormalCache* c, const Mesh* src, float group_eps) {
    if (!c || !src || src->nv <= 0 || src->nf <= 0) return -1;
    if (src->nv > MAXV || src->nf > MAXF) return -1;

    memset(c, 0, sizeof(*c));
    c->mesh = *src;
    c->group_eps = group_eps;
    c->position_groups = (group_eps >= 0.0f);

    Mesh* m = &c->mesh;
    for (int f = 0; f < m->nf; f++) {
        Tri t = m->f[f];
        if (t.a<0 || t.b<0 || t.c<0 || t.a>=m->nv || t.b>=m->nv || t.c>=m->nv)
            return -1;

        c->face_n[f] = cross(sub(m->v[t.b], m->v[t.a]),
                             sub(m->v[t.c], m->v[t.a]));
        c->accum[t.a] = add3(c->accum[t.a], c->face_n[f]);
        c->accum[t.b] = add3(c->accum[t.b], c->face_n[f]);
        c->accum[t.c] = add3(c->accum[t.c], c->face_n[f]);

        int ids[3] = { t.a, t.b, t.c };
        for (int k = 0; k < 3; k++) {
            int v = ids[k];
            int known = 0;
            for (int j = 0; j < c->incident_n[v]; j++) {
                if (c->incident[v][j] == f) { known = 1; break; }
            }
            /* Repeated corners give a zero face normal; one incidence is enough. */
            if (!known) c->incident[v][c->incident_n[v]++] = f;
        }
    }

    build_position_groups(m, group_eps, c->group, c->members);
    for (int g = 0; g < m->nv; g++) {
        if (c->members[g]) cache_recompute_group(c, g);
    }
    return 0;
}

static void cache_apply_face(NormalCache* c, int f, int add) {
    Tri t = c->mesh.f[f];
    V3 n = c->face_n[f];
    if (!add) n = scale3(n, -1.0f);
    c->accum[t.a] = add3(c->accum[t.a], n);
    c->accum[t.b] = add3(c->accum[t.b], n);
    c->accum[t.c] = add3(c->accum[t.c], n);
}

static unsigned group_closure(const NormalCache* c, unsigned seeds) {
    unsigned closure = 0;
    while (seeds) {
        int i = __builtin_ctz(seeds);
        closure |= c->members[c->group[i]];
        seeds &= seeds - 1u;
    }
    return closure;
}

/*
 * Incremental maintenance for exact co-location groups.  Exact equality is an
 * equivalence relation, so moving one vertex can only leave its old group and
 * optionally join one existing group.  Epsilon clustering is handled by a full
 * DSU rebuild because threshold chains can split when a bridge vertex moves.
 */
static void update_exact_position_group(NormalCache* c, int v) {
    int old_group = c->group[v];
    c->members[old_group] &= ~(1u << v);

    /* Keep the invariant that a non-empty group id is one of its members. */
    unsigned remaining = c->members[old_group];
    if (remaining && !(remaining & (1u << old_group))) {
        int new_group = __builtin_ctz(remaining);
        unsigned mask = remaining;
        while (mask) {
            int i = __builtin_ctz(mask);
            c->group[i] = new_group;
            mask &= mask - 1u;
        }
        c->members[new_group] = remaining;
        c->members[old_group] = 0;
    }

    int target = -1;
    for (int u = 0; u < c->mesh.nv; u++) {
        if (u == v) continue;
        if (c->mesh.v[u].x == c->mesh.v[v].x &&
            c->mesh.v[u].y == c->mesh.v[v].y &&
            c->mesh.v[u].z == c->mesh.v[v].z) {
            target = c->group[u];
            break;
        }
    }

    if (target >= 0) {
        c->group[v] = target;
        c->members[target] |= 1u << v;
    } else {
        c->group[v] = v;
        c->members[v] |= 1u << v;
    }
}

/*
 * Move one vertex and update only its dependency closure.
 *
 * Geometry changes face normals only on faces incident to v.  Every corner of
 * those faces is in the topological one-ring mask.  If normals are shared in
 * groups, a changed accumulator can affect every member of that corner's group;
 * regrouping can additionally revive the old group and expose the new group.
 * The returned mask is the exact set of output vertex normals rewritten.
 */
static unsigned normal_cache_move_vertex(NormalCache* c, int v,
                                         float x, float y, float z) {
    if (!c || v < 0 || v >= c->mesh.nv) return 0u;

    c->updated_faces = 0;
    unsigned changed = 1u << v;

    for (int k = 0; k < c->incident_n[v]; k++) {
        cache_apply_face(c, c->incident[v][k], 0 /* remove old */);
    }

    c->mesh.v[v] = (V3){x, y, z};

    for (int k = 0; k < c->incident_n[v]; k++) {
        int f = c->incident[v][k];
        Tri t = c->mesh.f[f];
        c->face_n[f] = cross(sub(c->mesh.v[t.b], c->mesh.v[t.a]),
                             sub(c->mesh.v[t.c], c->mesh.v[t.a]));
        cache_apply_face(c, f, 1 /* add new */);
        c->updated_faces++;
        changed |= (1u << t.a) | (1u << t.b) | (1u << t.c);
    }

    /* Expand through the old grouping before rebuilding it. */
    unsigned old_closure = group_closure(c, changed);

    if (c->position_groups) {
        if (c->group_eps == 0.0f)
            update_exact_position_group(c, v);
        else
            build_position_groups(&c->mesh, c->group_eps, c->group, c->members);
    }

    /* Old-group seeds can now belong to changed new groups. */
    unsigned dirty = group_closure(c, old_closure);

    unsigned recomputed_groups = 0;
    unsigned pending = dirty;
    while (pending) {
        int i = __builtin_ctz(pending);
        pending &= pending - 1u;

        int g = c->group[i];
        unsigned gbit = 1u << g;
        if (recomputed_groups & gbit) continue;
        recomputed_groups |= gbit;
        cache_recompute_group(c, g);
    }

    c->dirty_mask = dirty;
    return dirty;
}

/* Explicit normal-group change for callers that own smoothing/shared groups. */
static unsigned normal_cache_set_group(NormalCache* c, int v, int g) {
    if (!c || v < 0 || v >= c->mesh.nv || g < 0 || g >= c->mesh.nv) return 0u;

    int old_g = c->group[v];
    unsigned old_members = c->members[old_g];
    if (old_g == g) return old_members;

    c->members[old_g] &= ~(1u << v);
    c->group[v] = g;
    c->members[g] |= 1u << v;

    cache_recompute_group(c, old_g);
    cache_recompute_group(c, g);

    unsigned dirty = old_members | c->members[g];
    c->dirty_mask = dirty;
    return dirty;
}

/*
 * Deliberately wrong boundary, kept only for the regression counterexample.
 * It performs the same face/accumulator update but rewrites just the
 * topological one-ring indices and ignores normal-group aliases/regrouping.
 */
static unsigned counterexample_naive_one_ring_move(NormalCache* c, int v,
                                                   float x, float y, float z) {
    unsigned ring = 1u << v;

    for (int k = 0; k < c->incident_n[v]; k++) {
        cache_apply_face(c, c->incident[v][k], 0);
    }
    c->mesh.v[v] = (V3){x, y, z};
    for (int k = 0; k < c->incident_n[v]; k++) {
        int f = c->incident[v][k];
        Tri t = c->mesh.f[f];
        c->face_n[f] = cross(sub(c->mesh.v[t.b], c->mesh.v[t.a]),
                             sub(c->mesh.v[t.c], c->mesh.v[t.a]));
        cache_apply_face(c, f, 1);
        ring |= (1u << t.a) | (1u << t.b) | (1u << t.c);
    }

    for (int i = 0; i < c->mesh.nv; i++) {
        if (ring & (1u << i)) c->normal[i] = norm(c->accum[i]);
    }
    c->dirty_mask = ring;
    return ring;
}

static int v3_near(V3 a, V3 b, float eps) {
    return fabsf(a.x-b.x) <= eps && fabsf(a.y-b.y) <= eps && fabsf(a.z-b.z) <= eps;
}

static int cache_matches_grouped_full(const NormalCache* c, float eps) {
    int rebuilt_group[MAXV];
    unsigned rebuilt_members[MAXV];
    const int* group = c->group;

    /* Do not trust the cache's own grouping when checking a regrouping move.
     * Group ids are arbitrary labels; compare the induced equivalence relation. */
    if (c->position_groups) {
        build_position_groups(&c->mesh, c->group_eps,
                              rebuilt_group, rebuilt_members);
        group = rebuilt_group;
        for (int i = 0; i < c->mesh.nv; i++) {
            for (int j = i + 1; j < c->mesh.nv; j++) {
                int same_cache = c->group[i] == c->group[j];
                int same_full = rebuilt_group[i] == rebuilt_group[j];
                if (same_cache != same_full) return 0;
            }
        }
    }

    V3 full[MAXV];
    grouped_vertex_normals(&c->mesh, group, full);
    for (int i = 0; i < c->mesh.nv; i++) {
        if (!v3_near(c->normal[i], full[i], eps)) return 0;
    }
    return 1;
}

static int cache_matches_full(const NormalCache* c, float eps) {
    V3 full[MAXV];
    if (vertex_normals(&c->mesh, full) < 0) return 0;
    for (int i = 0; i < c->mesh.nv; i++) {
        if (!v3_near(c->normal[i], full[i], eps)) return 0;
    }
    return 1;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/*
 * Two disconnected reference triangles plus an unrelated triangle:
 *   0,1,2: moved vertex 0 starts at the origin
 *   3,4,5: vertex 3 is initially co-located/grouped with 0
 *   6,7,8: vertex 6 is far away at (2,0,0); 0 moves onto 6
 *   9..11: untouched control triangle
 *
 * Thus moving 0 splits old group {0,3} and merges new group {0,6}.
 */
static void make_group_spread_mesh(Mesh* m) {
    memset(m, 0, sizeof(*m));
    add_v(m, 0.0f, 0.0f, 0.0f);  /* 0, moved */
    add_v(m, 0.0f, 1.0f, 0.0f);  /* 1 */
    add_v(m, 0.0f, 0.0f, 1.0f);  /* 2 */

    add_v(m, 0.0f, 0.0f, 0.0f);  /* 3, old normal-group partner */
    add_v(m, 1.0f, 0.0f, 0.0f);  /* 4 */
    add_v(m, 1.0f, 1.0f, 0.0f);  /* 5 */

    add_v(m, 2.0f, 0.0f, 0.0f);  /* 6, distant new normal-group partner */
    add_v(m, 3.0f, 0.0f, 0.0f);  /* 7 */
    add_v(m, 2.0f, 1.0f, 0.0f);  /* 8 */

    add_v(m, 5.0f, 0.0f, 0.0f);  /* 9, untouched control */
    add_v(m, 6.0f, 0.0f, 0.0f);  /* 10 */
    add_v(m, 5.0f, 1.0f, 0.0f);  /* 11 */

    add_f(m, 0, 1, 2);
    add_f(m, 3, 4, 5);
    add_f(m, 6, 7, 8);
    add_f(m, 9, 10, 11);
}

int main(void) {
    Mesh m; memset(&m, 0, sizeof(m));
    add_v(&m, 0,0,0); add_v(&m, 1,0,0); add_v(&m, 0,1,0);
    add_f(&m, 0, 1, 2);
    V3 n[MAXV];
    int bad = vertex_normals(&m, n);
    check("no_degenerate", bad == 0);
    check("unit_length", fabsf(len(n[0]) - 1.0f) < 1e-5f);
    check("points_z", fabsf(fabsf(n[0].z) - 1.0f) < 1e-5f);

    /* two triangles of different area sharing an edge: weighting must favour big one */
    Mesh w; memset(&w, 0, sizeof(w));
    add_v(&w, 0,0,0); add_v(&w, 1,0,0); add_v(&w, 0,1,0);
    add_v(&w, 0,0,2);                      /* tall triangle */
    add_f(&w, 0, 1, 2);                    /* small, normal +z */
    add_f(&w, 0, 1, 3);                    /* large, normal -y-ish */
    V3 nw[MAXV];
    vertex_normals(&w, nw);
    V3 big = cross(sub(w.v[1], w.v[0]), sub(w.v[3], w.v[0]));
    V3 small = cross(sub(w.v[1], w.v[0]), sub(w.v[2], w.v[0]));
    check("bigger_dominates", len(big) > len(small));
    /* shared vertex 0 normal should tilt toward big triangle's normal */
    check("tilt_toward_big", fabsf(nw[0].y) > fabsf(nw[0].z) * 0.5f);

    Mesh empty; memset(&empty, 0, sizeof(empty));
    add_v(&empty, 0,0,0);
    check("no_faces_rejected", vertex_normals(&empty, n) == -1);

    /* Ungrouped incremental cache on a small grid: must equal full recompute. */
    Mesh grid; memset(&grid, 0, sizeof(grid));
    enum { GW = 4, GH = 4 };
    for (int y = 0; y < GH; y++)
        for (int x = 0; x < GW; x++)
            add_v(&grid, (float)x, (float)y, 0.0f);
    for (int y = 0; y < GH - 1; y++) {
        for (int x = 0; x < GW - 1; x++) {
            int a = y * GW + x;
            int b = a + 1;
            int d = a + GW;
            int e = d + 1;
            add_f(&grid, a, b, e);
            add_f(&grid, a, e, d);
        }
    }

    int grid_matches = 1;
    int face_count_ok = 1;
    int dirty_ring_ok = 1;
    for (int v = 0; v < grid.nv; v++) {
        NormalCache gc;
        if (normal_cache_init(&gc, &grid, -1.0f) != 0) {
            grid_matches = 0;
            break;
        }

        float px = gc.mesh.v[v].x + 0.25f * (float)((v % 3) - 1);
        float py = gc.mesh.v[v].y + 0.20f * (float)(((v + 1) % 3) - 1);
        float pz = 0.35f + 0.07f * (float)(v % 5);
        unsigned dirty = normal_cache_move_vertex(&gc, v, px, py, pz);

        if (!cache_matches_full(&gc, 2e-5f)) grid_matches = 0;

        int expected_faces = 0;
        unsigned expected_dirty = 1u << v;
        for (int f = 0; f < grid.nf; f++) {
            Tri t = grid.f[f];
            if (t.a == v || t.b == v || t.c == v) {
                expected_faces++;
                expected_dirty |= (1u << t.a) | (1u << t.b) | (1u << t.c);
            }
        }
        if (gc.updated_faces != expected_faces) face_count_ok = 0;
        if (dirty != expected_dirty) dirty_ring_ok = 0;
    }
    check("incremental_matches_full_grid", grid_matches);
    check("incremental_recomputes_only_incident_faces", face_count_ok);
    check("ungrouped_dirty_set_is_one_ring", dirty_ring_ok);

    /* Repeated moves on one cache must stay equal to repeated full recomputes. */
    NormalCache sequential;
    int sequential_ok = normal_cache_init(&sequential, &grid, -1.0f) == 0;
    if (sequential_ok) {
        for (int step = 0; step < 24 && sequential_ok; step++) {
            int v = step % grid.nv;
            float z = 0.15f + 0.05f * (float)(step % 7);
            normal_cache_move_vertex(&sequential, v,
                                     sequential.mesh.v[v].x,
                                     sequential.mesh.v[v].y, z);
            sequential_ok = cache_matches_full(&sequential, 1e-4f);
        }
    }
    check("repeated_incremental_moves_match_full", sequential_ok);

    /* Counterexample: a topological one-ring update misses group spread. */
    Mesh spread;
    make_group_spread_mesh(&spread);

    NormalCache baseline_groups;
    int spread_init_ok = normal_cache_init(&baseline_groups, &spread, 0.0f) == 0;
    V3 before[MAXV];
    memcpy(before, baseline_groups.normal, sizeof(before));

    NormalCache ring_only;
    int ring_init_ok = normal_cache_init(&ring_only, &spread, 0.0f) == 0;
    unsigned ring_mask = counterexample_naive_one_ring_move(&ring_only, 0, 2,0,0);

    unsigned dirty = 0;
    NormalCache incremental_groups;
    int inc_init_ok = normal_cache_init(&incremental_groups, &spread, 0.0f) == 0;
    if (inc_init_ok) dirty = normal_cache_move_vertex(&incremental_groups, 0, 2,0,0);

    V3 expected_grouped[MAXV];
    int expected_group[MAXV];
    unsigned expected_members[MAXV];
    build_position_groups(&incremental_groups.mesh, 0.0f,
                          expected_group, expected_members);
    grouped_vertex_normals(&incremental_groups.mesh,
                           expected_group, expected_grouped);

    int ring_excludes_far = (ring_mask & (1u << 6)) == 0;
    int far_true_changed = !v3_near(before[6], expected_grouped[6], 1e-4f);
    int ring_left_far_stale = v3_near(ring_only.normal[6], before[6], 1e-6f);
    int ring_wrong_vs_full = !v3_near(ring_only.normal[6], expected_grouped[6], 1e-4f);
    check("spread_fixture_init", spread_init_ok && ring_init_ok && inc_init_ok);
    check("counterexample_ring_excludes_far", ring_excludes_far);
    check("counterexample_far_normal_changes", far_true_changed);
    check("counterexample_one_ring_misses_spread",
          ring_excludes_far && far_true_changed &&
          ring_left_far_stale && ring_wrong_vs_full);

    /* The correct boundary catches both old-group split and new-group merge. */
    unsigned expected_group_dirty =
        (1u << 0) | (1u << 1) | (1u << 2) | (1u << 3) | (1u << 6);
    check("grouped_dirty_closure_exact", dirty == expected_group_dirty);
    check("incremental_grouped_matches_full",
          cache_matches_grouped_full(&incremental_groups, 2e-5f));
    check("incremental_updates_old_group_partner",
          (dirty & (1u << 3)) &&
          !v3_near(incremental_groups.normal[3], before[3], 1e-4f));
    check("incremental_updates_distant_group_partner",
          (dirty & (1u << 6)) &&
          !v3_near(incremental_groups.normal[6], before[6], 1e-4f));
    check("group_spread_spares_unrelated",
          (dirty & ((1u << 9) | (1u << 10) | (1u << 11))) == 0);
    check("grouped_move_recomputes_one_face",
          incremental_groups.updated_faces == 1);

    /* The same closure works for arbitrary explicit groups, not just co-location. */
    NormalCache explicit_alias;
    int explicit_alias_ok = normal_cache_init(&explicit_alias, &spread, -1.0f) == 0;
    unsigned explicit_group_dirty = normal_cache_set_group(&explicit_alias, 9, 1);
    check("explicit_group_merge_marks_members",
          explicit_alias_ok &&
          (explicit_group_dirty & ((1u << 1) | (1u << 9))) ==
              ((1u << 1) | (1u << 9)));

    unsigned explicit_move_dirty =
        normal_cache_move_vertex(&explicit_alias, 0, 2.0f, 0.0f, 0.0f);
    check("explicit_group_alias_spreads_beyond_ring",
          (explicit_move_dirty & (1u << 9)) &&
          !(explicit_move_dirty & (1u << 6)) &&
          cache_matches_grouped_full(&explicit_alias, 2e-5f));

    /* Repeated split/merge moves must keep the incremental groups canonical. */
    NormalCache grouped_sequence;
    int grouped_sequence_ok =
        normal_cache_init(&grouped_sequence, &spread, 0.0f) == 0;
    if (grouped_sequence_ok) {
        const float path[][3] = {
            {2,0,0}, {0,0,0}, {2,0,0}, {4,0,0}, {0,0,0}
        };
        for (int step = 0; step < 5; step++) {
            normal_cache_move_vertex(&grouped_sequence, 0,
                                     path[step][0], path[step][1], path[step][2]);
            if (!cache_matches_grouped_full(&grouped_sequence, 1e-4f)) {
                grouped_sequence_ok = 0;
                break;
            }
        }
    }
    check("repeated_group_moves_match_full", grouped_sequence_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
