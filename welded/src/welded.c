
/* welded: order-independent vertex welding via spatial hash grid + union-find,
   with face remap and degenerate face removal. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXV 100000
#define MAXF 100000

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct {
    V3 v[MAXV]; int nv;
    Tri f[MAXF]; int nf;
    int map[MAXV];      /* old index -> welded index */
    int welded_count;
} Mesh;

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); for (int i = 0; i < MAXV; i++) m->map[i] = -1; }

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

/* ---- union-find: equivalence classes are connected components of the
   "distance <= eps" graph, which depends only on the point set, not order ---- */

static int dsu_find(int* parent, int x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}

static void dsu_union(int* parent, int* size, int a, int b) {
    a = dsu_find(parent, a);
    b = dsu_find(parent, b);
    if (a == b) return;
    if (size[a] < size[b]) { int t = a; a = b; b = t; }
    parent[b] = a;
    size[a] += size[b];
}

/* lexicographic order on stored float bits: deterministic, order-free */
static int v3_lex(const V3* a, const V3* b) {
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    if (a->z < b->z) return -1;
    if (a->z > b->z) return 1;
    return 0;
}

/* spatial hash: cells of side eps, chained vertices per occupied cell */
static unsigned cell_hash(int cx, int cy, int cz) {
    return ((unsigned)cx * 73856093u)
         ^ ((unsigned)cy * 19349663u)
         ^ ((unsigned)cz * 83492791u);
}

static const V3* g_canon_pos;   /* qsort context */
static int* g_canon_min;        /* root -> lexicographically smallest member */
static int canon_cmp(const void* pa, const void* pb) {
    const V3* a = &g_canon_pos[g_canon_min[*(const int*)pa]];
    const V3* b = &g_canon_pos[g_canon_min[*(const int*)pb]];
    return v3_lex(a, b);
}

/* weld vertices closer than eps, remap faces, drop degenerate.
   The partition is independent of vertex input order: two vertices end up in
   the same class iff connected by a chain of <=eps edges. Returns class count. */
static int mesh_weld(Mesh* m, float eps) {
    if (eps <= 0.0f) return -1;
    int n = m->nv;
    if (n == 0) { m->welded_count = 0; return 0; }

    int* parent = (int*)malloc(sizeof(int) * (size_t)n);
    int* size   = (int*)malloc(sizeof(int) * (size_t)n);
    int* next   = (int*)malloc(sizeof(int) * (size_t)n);
    int* cmin   = (int*)malloc(sizeof(int) * (size_t)n);
    if (!parent || !size || !next || !cmin) {
        free(parent); free(size); free(next); free(cmin);
        return -1;
    }

    /* table sized for load factor <= 0.25 */
    unsigned cap = 16;
    while (cap < (unsigned)n * 4u) cap <<= 1;
    int* head = (int*)malloc(sizeof(int) * cap);
    if (!head) { free(parent); free(size); free(next); free(cmin); return -1; }
    memset(head, -1, sizeof(int) * cap);
    unsigned mask = cap - 1;

    for (int i = 0; i < n; i++) { parent[i] = i; size[i] = 1; }

    float eps2 = eps * eps;
    for (int i = 0; i < n; i++) {
        int cx = (int)floorf(m->v[i].x / eps);
        int cy = (int)floorf(m->v[i].y / eps);
        int cz = (int)floorf(m->v[i].z / eps);

        /* any point within eps differs by at most one cell index per axis,
           so the 27 neighboring cells exhaust every possible <=eps pair */
        for (int dz = -1; dz <= 1; dz++)
        for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++) {
            unsigned h = cell_hash(cx + dx, cy + dy, cz + dz) & mask;
            for (int j = head[h]; j >= 0; j = next[j]) {
                float ex = m->v[i].x - m->v[j].x;
                float ey = m->v[i].y - m->v[j].y;
                float ez = m->v[i].z - m->v[j].z;
                if (ex*ex + ey*ey + ez*ez <= eps2)
                    dsu_union(parent, size, i, j);
            }
        }

        unsigned h = cell_hash(cx, cy, cz) & mask;
        next[i] = head[h];
        head[h] = i;
    }

    /* canonical label per component: order roots by the lexicographically
       smallest member position. Distinct components cannot share a position
       (identical positions are distance 0 and would have been unioned), so
       this ordering has no order-dependent ties. */
    for (int i = 0; i < n; i++) cmin[i] = i;
    for (int i = 0; i < n; i++) {
        int r = dsu_find(parent, i);
        if (v3_lex(&m->v[i], &m->v[cmin[r]]) < 0) cmin[r] = i;
    }
    int k = 0;
    for (int i = 0; i < n; i++) if (dsu_find(parent, i) == i) k++;
    int* roots = (int*)malloc(sizeof(int) * (size_t)k);
    int* label = (int*)malloc(sizeof(int) * (size_t)n);
    V3* uniq  = (V3*)malloc(sizeof(V3) * (size_t)k);
    if (!roots || !label || !uniq) {
        free(parent); free(size); free(next); free(cmin); free(head);
        free(roots); free(label); free(uniq);
        return -1;
    }
    k = 0;
    for (int i = 0; i < n; i++) if (dsu_find(parent, i) == i) roots[k++] = i;
    g_canon_pos = m->v;
    g_canon_min = cmin;
    qsort(roots, (size_t)k, sizeof(int), canon_cmp);

    for (int t = 0; t < k; t++) {
        int r = roots[t];
        label[r] = t;
        uniq[t] = m->v[cmin[r]];      /* representative: canonical member */
    }
    for (int i = 0; i < n; i++) m->map[i] = label[dsu_find(parent, i)];

    int nf = 0;
    for (int i = 0; i < m->nf; i++) {
        int a = m->map[m->f[i].a], b = m->map[m->f[i].b], c = m->map[m->f[i].c];
        if (a == b || b == c || a == c) continue;
        m->f[nf].a = a; m->f[nf].b = b; m->f[nf].c = c; nf++;
    }
    m->nf = nf;
    memcpy(m->v, uniq, sizeof(V3) * (size_t)k);
    m->nv = k;
    m->welded_count = k;

    free(parent); free(size); free(next); free(cmin); free(head);
    free(roots); free(label); free(uniq);
    return k;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---------------- order-invariance test machinery ---------------- */

static unsigned rng_state;
static unsigned rng_next(void) {
    unsigned x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_state = x;
    return x;
}
static void rng_seed(unsigned s) { rng_state = s ? s : 1; }
static float rng_coord(void) { return (float)(rng_next() % 10000) / 1000.0f; } /* [0,10) */
static void shuffle_perm(int* p, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = (int)(rng_next() % (unsigned)(i + 1));
        int t = p[i]; p[i] = p[j]; p[j] = t;
    }
}

/* build a mesh from pts/faces in permuted slot order, weld, and report the
   welded label of every ORIGINAL vertex id (order-free comparison anchor) */
static int run_welded(Mesh* m, const V3* pts, int n, const int* perm,
                      const Tri* faces, int nface, float eps, int* label_orig) {
    int* inv = (int*)malloc(sizeof(int) * (size_t)n);
    if (!inv) return -1;
    mesh_init(m);
    for (int s = 0; s < n; s++) {
        inv[perm[s]] = s;
        mesh_add_vertex(m, pts[perm[s]].x, pts[perm[s]].y, pts[perm[s]].z);
    }
    for (int f = 0; f < nface; f++)
        mesh_add_tri(m, inv[faces[f].a], inv[faces[f].b], inv[faces[f].c]);
    int nu = mesh_weld(m, eps);
    for (int o = 0; o < n; o++) label_orig[o] = m->map[inv[o]];
    free(inv);
    return nu;
}

/* reference: brute-force O(n^2) enumeration of every <=eps pair, plain DSU */
static void ref_partition(const V3* pts, int n, float eps, int* root_out) {
    int* par = (int*)malloc(sizeof(int) * (size_t)n);
    for (int i = 0; i < n; i++) par[i] = i;
    float eps2 = eps * eps;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            float dx = pts[i].x - pts[j].x;
            float dy = pts[i].y - pts[j].y;
            float dz = pts[i].z - pts[j].z;
            if (dx*dx + dy*dy + dz*dz <= eps2) {
                int a = i, b = j;
                while (par[a] != a) a = par[a];
                while (par[b] != b) b = par[b];
                if (a != b) par[b] = a;
            }
        }
    for (int i = 0; i < n; i++) {
        int r = i; while (par[r] != r) r = par[r];
        root_out[i] = r;
    }
    free(par);
}

int main(void) {
    /* ---- original behavior tests ---- */
    static Mesh m;
    mesh_init(&m);
    mesh_add_vertex(&m, 0, 0, 0);
    mesh_add_vertex(&m, 1, 0, 0);
    mesh_add_vertex(&m, 0, 1, 0);
    mesh_add_vertex(&m, 0, 0, 0);        /* duplicate of 0 */
    mesh_add_vertex(&m, 1.0001f, 0, 0);  /* near duplicate of 1 */
    check("added", m.nv == 5);

    check("tri_ok", mesh_add_tri(&m, 0, 1, 2) == 0);
    check("tri_dup_rejected", mesh_add_tri(&m, 0, 0, 2) == -1);
    check("tri_range_rejected", mesh_add_tri(&m, 0, 1, 99) == -1);
    mesh_add_tri(&m, 3, 4, 2);           /* welds to same as first: survives */
    mesh_add_tri(&m, 0, 3, 1);           /* 0 and 3 weld together: degenerate, dropped */

    int nu = mesh_weld(&m, 1e-3f);
    check("weld_count", nu == 3);
    check("map_zero", m.map[0] == m.map[3]);
    check("map_one", m.map[1] == m.map[4]);
    check("degenerate_dropped", m.nf == 2);

    check("bad_eps", mesh_weld(&m, 0.0f) == -1 && mesh_weld(&m, -1.0f) == -1);

    static Mesh e;
    mesh_init(&e);
    mesh_add_vertex(&e, 0,0,0); mesh_add_vertex(&e, 1,0,0); mesh_add_vertex(&e, 0,1,0);
    mesh_add_tri(&e, 0, 1, 2);
    check("no_weld_needed", mesh_weld(&e, 1e-4f) == 3 && e.nf == 1);

    static Mesh empty_mesh;
    mesh_init(&empty_mesh);
    check("empty_weld", mesh_weld(&empty_mesh, 1e-3f) == 0);

    /* ---- chain case: greedy first-fit would split this differently by order;
       connected-component welding must give one class in both orders ---- */
    V3 chain[3] = {{0,0,0},{0.8f,0,0},{1.6f,0,0}}; /* consecutive dist 0.8 <= 1 */
    int fwd[3] = {0,1,2}, rev[3] = {2,1,0};
    static Mesh cm;
    int lab_fwd[3], lab_rev[3];
    int nf_fwd = run_welded(&cm, chain, 3, fwd, NULL, 0, 1.0f, lab_fwd);
    int nf_rev = run_welded(&cm, chain, 3, rev, NULL, 0, 1.0f, lab_rev);
    check("chain_single_class", nf_fwd == 1 && nf_rev == 1);
    check("chain_order_invariant",
          lab_fwd[0] == lab_fwd[1] && lab_fwd[1] == lab_fwd[2] &&
          lab_rev[0] == lab_rev[1] && lab_rev[1] == lab_rev[2]);

    /* ---- shuffled-input invariance on random data with duplicates + faces ---- */
    enum { N = 300, NFACE = 40, RUNS = 5 };
    static V3 pts[N];
    static Tri faces[NFACE];
    rng_seed(0xC0FFEEu);
    for (int i = 0; i < N; i++) pts[i] = (V3){ rng_coord(), rng_coord(), rng_coord() };
    pts[7]  = pts[0];                 /* exact duplicates */
    pts[88] = pts[0];
    pts[201] = pts[42];
    for (int f = 0; f < NFACE; f++) {
        int a = (int)(rng_next() % N);
        int b = (int)(rng_next() % N); while (b == a) b = (int)(rng_next() % N);
        int c = (int)(rng_next() % N); while (c == a || c == b) c = (int)(rng_next() % N);
        faces[f] = (Tri){ a, b, c };
    }

    static int perm[N], base_label[N], cur_label[N];
    static Tri base_faces[NFACE];
    static V3 base_pos[N];
    int base_nu = 0, base_nf = 0, base_nv = 0;
    int shuffle_ok = 1, faces_ok = 1, positions_ok = 1;
    for (int run = 0; run < RUNS; run++) {
        for (int i = 0; i < N; i++) perm[i] = i;
        if (run > 0) shuffle_perm(perm, N);
        int nu2 = run_welded(&cm, pts, N, perm, faces, NFACE, 0.25f, cur_label);
        if (run == 0) {
            base_nu = nu2; base_nf = cm.nf; base_nv = cm.nv;
            memcpy(base_label, cur_label, sizeof(base_label));
            memcpy(base_faces, cm.f, sizeof(Tri) * (size_t)cm.nf);
            memcpy(base_pos, cm.v, sizeof(V3) * (size_t)cm.nv);
        } else {
            if (nu2 != base_nu) shuffle_ok = 0;
            for (int i = 0; i < N; i++)
                if (cur_label[i] != base_label[i]) shuffle_ok = 0;
            if (cm.nf != base_nf) faces_ok = 0;
            for (int i = 0; i < base_nf; i++)
                if (memcmp(&cm.f[i], &base_faces[i], sizeof(Tri)) != 0) faces_ok = 0;
            if (cm.nv != base_nv) positions_ok = 0;
            for (int i = 0; i < base_nv; i++)
                if (memcmp(&cm.v[i], &base_pos[i], sizeof(V3)) != 0) positions_ok = 0;
        }
    }
    check("shuffle_partition_invariant", shuffle_ok);
    check("shuffle_face_invariant", faces_ok);
    check("shuffle_position_invariant", positions_ok);

    /* ---- brute-force O(n^2) reference: bucket grid must miss no pair ---- */
    enum { RN = 60 };
    static V3 rpts[RN];
    rng_seed(0xBEEFu);
    for (int i = 0; i < RN; i++) rpts[i] = (V3){ rng_coord(), rng_coord(), rng_coord() };
    rpts[33] = rpts[5]; rpts[59] = rpts[12];
    static int ref_root[RN], bucket_lab[RN];
    ref_partition(rpts, RN, 0.25f, ref_root);
    static Mesh rm;
    for (int i = 0; i < RN; i++) perm[i] = i;
    rng_seed(7u); shuffle_perm(perm, RN);
    run_welded(&rm, rpts, RN, perm, NULL, 0, 0.25f, bucket_lab);
    int ref_ok = 1;
    for (int i = 0; i < RN && ref_ok; i++)
        for (int j = i + 1; j < RN; j++) {
            int same_ref = ref_root[i] == ref_root[j];
            int same_bucket = bucket_lab[i] == bucket_lab[j];
            if (same_ref != same_bucket) { ref_ok = 0; break; }
        }
    check("matches_bruteforce_pairs", ref_ok);

    /* ---- scale: timing at 10k and 20k should grow ~linearly, not quadratically ---- */
    static V3 big[20000];
    static Mesh bm;
    rng_seed(0x5151u);
    for (int i = 0; i < 20000; i++) big[i] = (V3){ rng_coord(), rng_coord(), rng_coord() };
    double ms10k = 0, ms20k = 0; int w10k = 0, w20k = 0;
    for (int size = 10000; size <= 20000; size += 10000) {
        mesh_init(&bm);
        for (int i = 0; i < size; i++)
            mesh_add_vertex(&bm, big[i].x, big[i].y, big[i].z);
        clock_t t0 = clock();
        int w = mesh_weld(&bm, 0.05f);
        double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
        if (size == 10000) { ms10k = ms; w10k = w; } else { ms20k = ms; w20k = w; }
    }
    check("scale_classes_sane", w10k > 0 && w10k <= 10000 && w20k > 0 && w20k <= 20000);
    printf("weld_10k_ms=%.3f welded_10k=%d\n", ms10k, w10k);
    printf("weld_20k_ms=%.3f welded_20k=%d\n", ms20k, w20k);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
