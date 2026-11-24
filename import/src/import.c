
/* import: OBJ model import pipeline, split into timed stages.
 *
 *   read  -> parse -> weld -> index
 *
 * The old loader welded vertices *inside* the parse loop, so a single
 * timer around "load" attributed the O(n^2) weld to parsing.  This
 * module keeps the stages separate and times each one.
 *
 *   ./import file.obj          split timing, naive O(n^2) weld
 *   ./import --fast file.obj   split timing, hash-grid weld
 *   ./import --legacy file.obj old loader: one timer, weld inlined in parse
 *   ./import --selftest        correctness checks (naive vs grid weld)
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct { V3* v; int nv, capv; Tri* f; int nf, capf; } Mesh;

/* vertex->face adjacency, CSR */
typedef struct { int* off; int* adj; } Index;

typedef struct {
    double read_ms, parse_ms, weld_ms, index_ms;
} Timings;

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e3 + ts.tv_nsec / 1e6;
}

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

static void mesh_free(Mesh* m) { free(m->v); free(m->f); mesh_init(m); }

static int mesh_push_v(Mesh* m, float x, float y, float z) {
    if (m->nv == m->capv) {
        int cap = m->capv ? m->capv * 2 : 1024;
        V3* p = realloc(m->v, sizeof(V3) * cap);
        if (!p) return -1;
        m->v = p; m->capv = cap;
    }
    m->v[m->nv].x = x; m->v[m->nv].y = y; m->v[m->nv].z = z;
    return m->nv++;
}

static int mesh_push_f(Mesh* m, int a, int b, int c) {
    if (m->nf == m->capf) {
        int cap = m->capf ? m->capf * 2 : 1024;
        Tri* p = realloc(m->f, sizeof(Tri) * cap);
        if (!p) return -1;
        m->f = p; m->capf = cap;
    }
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* ---- stage 1: read ---------------------------------------------------- */

static char* read_file(const char* path, long* out_size) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return NULL;
    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char* buf = malloc(size + 1);
    if (!buf) { fclose(fp); return NULL; }
    if (fread(buf, 1, size, fp) != (size_t)size) { free(buf); fclose(fp); return NULL; }
    buf[size] = '\0';
    fclose(fp);
    *out_size = size;
    return buf;
}

/* ---- stage 2: parse (pure tokenizing, no welding) --------------------- */

static int parse_obj(Mesh* m, char* buf) {
    char* p = buf;
    while (*p) {
        if (p[0] == 'v' && p[1] == ' ') {
            char* end;
            float x = strtof(p + 2, &end);
            float y = strtof(end, &end);
            float z = strtof(end, &end);
            if (mesh_push_v(m, x, y, z) < 0) return -1;
        } else if (p[0] == 'f' && p[1] == ' ') {
            char* end;
            /* strtol stops at '/' in "a/b/c" forms; vertex index is all we need */
            long a = strtol(p + 2, &end, 10);
            long b = strtol(end, &end, 10);
            long c = strtol(end, &end, 10);
            if (mesh_push_f(m, (int)a - 1, (int)b - 1, (int)c - 1) < 0) return -1;
        }
        while (*p && *p != '\n') p++;
        if (*p) p++;
    }
    return 0;
}

/* ---- stage 3a: naive weld, same algorithm as welded.c ----------------- */

/* remap faces through map, drop degenerate, compact vertices to uniques */
static void weld_finish(Mesh* m, const V3* uniq, int nu, const int* map) {
    int nf = 0;
    for (int i = 0; i < m->nf; i++) {
        int a = map[m->f[i].a], b = map[m->f[i].b], c = map[m->f[i].c];
        if (a == b || b == c || a == c) continue;
        m->f[nf].a = a; m->f[nf].b = b; m->f[nf].c = c; nf++;
    }
    m->nf = nf;
    memcpy(m->v, uniq, sizeof(V3) * nu);
    m->nv = nu;
}

static int weld_naive(Mesh* m, float eps) {
    V3* uniq = malloc(sizeof(V3) * m->nv);
    int* map = malloc(sizeof(int) * m->nv);
    if (!uniq || !map) { free(uniq); free(map); return -1; }
    int nu = 0;
    for (int i = 0; i < m->nv; i++) {
        int found = -1;
        for (int j = 0; j < nu; j++) {
            float dx = m->v[i].x - uniq[j].x;
            float dy = m->v[i].y - uniq[j].y;
            float dz = m->v[i].z - uniq[j].z;
            if (dx*dx + dy*dy + dz*dz <= eps*eps) { found = j; break; }
        }
        if (found < 0) { uniq[nu] = m->v[i]; found = nu++; }
        map[i] = found;
    }
    weld_finish(m, uniq, nu, map);
    free(uniq); free(map);
    return nu;
}

/* ---- stage 3b: hash-grid weld, O(n) expected -------------------------- */

static unsigned cell_hash(int ix, int iy, int iz) {
    return (unsigned)(ix * 73856093) ^ (unsigned)(iy * 19349663) ^ (unsigned)(iz * 83492791);
}

static int weld_grid(Mesh* m, float eps) {
    V3* uniq = malloc(sizeof(V3) * m->nv);
    int* map = malloc(sizeof(int) * m->nv);
    int* next = malloc(sizeof(int) * m->nv);   /* chain per unique vertex */
    int size = 1024;
    while (size < m->nv * 2) size *= 2;
    int* head = malloc(sizeof(int) * size);
    if (!uniq || !map || !next || !head) { free(uniq); free(map); free(next); free(head); return -1; }
    for (int i = 0; i < size; i++) head[i] = -1;

    float inv = 1.0f / eps;
    int nu = 0;
    for (int i = 0; i < m->nv; i++) {
        int ix = (int)floorf(m->v[i].x * inv);
        int iy = (int)floorf(m->v[i].y * inv);
        int iz = (int)floorf(m->v[i].z * inv);
        int found = -1;
        /* any point within eps lies in the 27 neighboring cells */
        for (int dx = -1; dx <= 1 && found < 0; dx++)
        for (int dy = -1; dy <= 1 && found < 0; dy++)
        for (int dz = -1; dz <= 1 && found < 0; dz++) {
            unsigned h = cell_hash(ix + dx, iy + dy, iz + dz) & (size - 1);
            for (int j = head[h]; j >= 0; j = next[j]) {
                float ddx = m->v[i].x - uniq[j].x;
                float ddy = m->v[i].y - uniq[j].y;
                float ddz = m->v[i].z - uniq[j].z;
                if (ddx*ddx + ddy*ddy + ddz*ddz <= eps*eps) { found = j; break; }
            }
        }
        if (found < 0) {
            uniq[nu] = m->v[i];
            unsigned h = cell_hash(ix, iy, iz) & (size - 1);
            next[nu] = head[h]; head[h] = nu;
            found = nu++;
        }
        map[i] = found;
    }
    weld_finish(m, uniq, nu, map);
    free(uniq); free(map); free(next); free(head);
    return nu;
}

/* ---- stage 4: vertex->face adjacency index (CSR) ---------------------- */

static int build_index(const Mesh* m, Index* ix) {
    ix->off = calloc(m->nv + 1, sizeof(int));
    if (!ix->off) return -1;
    for (int i = 0; i < m->nf; i++) {
        ix->off[m->f[i].a + 1]++;
        ix->off[m->f[i].b + 1]++;
        ix->off[m->f[i].c + 1]++;
    }
    for (int i = 0; i < m->nv; i++) ix->off[i + 1] += ix->off[i];
    ix->adj = malloc(sizeof(int) * (ix->off[m->nv] ? ix->off[m->nv] : 1));
    if (!ix->adj) { free(ix->off); return -1; }
    int* cur = malloc(sizeof(int) * (m->nv ? m->nv : 1));
    if (!cur) { free(ix->off); free(ix->adj); return -1; }
    memcpy(cur, ix->off, sizeof(int) * m->nv);
    for (int i = 0; i < m->nf; i++) {
        ix->adj[cur[m->f[i].a]++] = i;
        ix->adj[cur[m->f[i].b]++] = i;
        ix->adj[cur[m->f[i].c]++] = i;
    }
    free(cur);
    return 0;
}

/* ---- legacy loader: weld inlined in the parse loop -------------------- */

/* What the old code did: every parsed vertex was checked against all
 * uniques so far, inside the same function (and the same timer) as
 * tokenizing.  Kept here to reproduce the original measurement. */
static int parse_obj_inline_weld(Mesh* m, char* buf, float eps) {
    V3* uniq = NULL; int nu = 0, capu = 0;
    int* map = NULL; int capm = 0;
    char* p = buf;
    while (*p) {
        if (p[0] == 'v' && p[1] == ' ') {
            char* end;
            float x = strtof(p + 2, &end);
            float y = strtof(end, &end);
            float z = strtof(end, &end);
            if (nu == capu) {
                capu = capu ? capu * 2 : 1024;
                uniq = realloc(uniq, sizeof(V3) * capu);
            }
            if (m->nv == capm) {
                capm = capm ? capm * 2 : 1024;
                map = realloc(map, sizeof(int) * capm);
            }
            int found = -1;
            for (int j = 0; j < nu; j++) {           /* O(n^2), hidden in "parse" */
                float dx = x - uniq[j].x, dy = y - uniq[j].y, dz = z - uniq[j].z;
                if (dx*dx + dy*dy + dz*dz <= eps*eps) { found = j; break; }
            }
            if (found < 0) { uniq[nu].x = x; uniq[nu].y = y; uniq[nu].z = z; found = nu++; }
            map[m->nv] = found;
            mesh_push_v(m, x, y, z);
        } else if (p[0] == 'f' && p[1] == ' ') {
            char* end;
            long a = strtol(p + 2, &end, 10);
            long b = strtol(end, &end, 10);
            long c = strtol(end, &end, 10);
            mesh_push_f(m, (int)a - 1, (int)b - 1, (int)c - 1);
        }
        while (*p && *p != '\n') p++;
        if (*p) p++;
    }
    weld_finish(m, uniq, nu, map);
    free(uniq); free(map);
    return nu;
}

/* ---- drivers ---------------------------------------------------------- */

static int run_split(const char* path, float eps, int fast) {
    Timings t;
    long size = 0;
    Mesh m; mesh_init(&m);
    Index ix;

    double t0 = now_ms();
    char* buf = read_file(path, &size);
    if (!buf) { fprintf(stderr, "read failed: %s\n", path); return 1; }
    t.read_ms = now_ms() - t0;

    t0 = now_ms();
    if (parse_obj(&m, buf) < 0) { fprintf(stderr, "parse failed\n"); return 1; }
    t.parse_ms = now_ms() - t0;
    free(buf);

    t0 = now_ms();
    int nu = fast ? weld_grid(&m, eps) : weld_naive(&m, eps);
    if (nu < 0) { fprintf(stderr, "weld failed\n"); return 1; }
    t.weld_ms = now_ms() - t0;

    t0 = now_ms();
    if (build_index(&m, &ix) < 0) { fprintf(stderr, "index failed\n"); return 1; }
    t.index_ms = now_ms() - t0;

    printf("mode=%s file=%s size=%.1fMB\n", fast ? "fast" : "naive", path, size / 1e6);
    printf("read_ms=%.1f\nparse_ms=%.1f\nweld_ms=%.1f\nindex_ms=%.1f\n",
           t.read_ms, t.parse_ms, t.weld_ms, t.index_ms);
    printf("total_ms=%.1f\n", t.read_ms + t.parse_ms + t.weld_ms + t.index_ms);
    printf("verts=%d faces=%d welded_to=%d\n", m.nv, m.nf, nu);

    free(ix.off); free(ix.adj);
    mesh_free(&m);
    return 0;
}

static int run_legacy(const char* path, float eps) {
    long size = 0;
    Mesh m; mesh_init(&m);
    Index ix;

    double t0 = now_ms();
    char* buf = read_file(path, &size);
    if (!buf) { fprintf(stderr, "read failed: %s\n", path); return 1; }
    double read_ms = now_ms() - t0;

    /* the old "parse" timer: tokenizing + inline weld, one number */
    t0 = now_ms();
    int nu = parse_obj_inline_weld(&m, buf, eps);
    double parse_ms = now_ms() - t0;
    free(buf);

    t0 = now_ms();
    if (build_index(&m, &ix) < 0) { fprintf(stderr, "index failed\n"); return 1; }
    double index_ms = now_ms() - t0;

    printf("mode=legacy file=%s size=%.1fMB\n", path, size / 1e6);
    printf("read_ms=%.1f\nparse_ms=%.1f\nindex_ms=%.1f\n", read_ms, parse_ms, index_ms);
    printf("total_ms=%.1f\n", read_ms + parse_ms + index_ms);
    printf("verts=%d faces=%d welded_to=%d\n", m.nv, m.nf, nu);

    free(ix.off); free(ix.adj);
    mesh_free(&m);
    return 0;
}

/* ---- selftest --------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void build_dup_mesh(Mesh* m) {
    mesh_init(m);
    mesh_push_v(m, 0, 0, 0);
    mesh_push_v(m, 1, 0, 0);
    mesh_push_v(m, 0, 1, 0);
    mesh_push_v(m, 0, 0, 0);            /* duplicate of 0 */
    mesh_push_v(m, 1.0001f, 0, 0);      /* near duplicate of 1 */
    mesh_push_f(m, 0, 1, 2);
    mesh_push_f(m, 3, 4, 2);            /* welds to same as first: survives */
    mesh_push_f(m, 0, 3, 1);            /* degenerate after weld: dropped */
}

static int selftest(void) {
    Mesh a, b;
    build_dup_mesh(&a);
    build_dup_mesh(&b);

    int na = weld_naive(&a, 1e-3f);
    check("naive_weld_count", na == 3);
    check("naive_degenerate_dropped", a.nf == 2);

    int nb = weld_grid(&b, 1e-3f);
    check("grid_weld_count", nb == na);
    check("grid_same_faces", b.nf == a.nf);
    check("grid_same_verts", b.nv == a.nv);

    /* grid vs naive on a denser deterministic mesh */
    Mesh c, d; mesh_init(&c); mesh_init(&d);
    for (int i = 0; i < 50; i++)
        for (int j = 0; j < 50; j++) {
            float x = i * 0.5f, y = j * 0.5f;
            mesh_push_v(&c, x, y, 0); mesh_push_v(&d, x, y, 0);
            mesh_push_v(&c, x, y, 0); mesh_push_v(&d, x, y, 0);  /* exact dup */
        }
    int nc = weld_naive(&c, 1e-3f);
    int nd = weld_grid(&d, 1e-3f);
    check("grid_matches_naive_dense", nc == nd && nc == 2500);

    Index ix;
    check("index_build", build_index(&a, &ix) == 0);
    long sum = ix.off[a.nv];
    check("index_adj_sum", sum == 3L * a.nf);
    check("index_first_face", ix.adj[ix.off[0]] == 0);
    free(ix.off); free(ix.adj);

    mesh_free(&a); mesh_free(&b); mesh_free(&c); mesh_free(&d);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    float eps = 1e-3f;
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc == 3 && strcmp(argv[1], "--fast") == 0) return run_split(argv[2], eps, 1);
    if (argc == 3 && strcmp(argv[1], "--legacy") == 0) return run_legacy(argv[2], eps);
    if (argc == 2) return run_split(argv[1], eps, 0);
    fprintf(stderr, "usage: import [--fast|--legacy] file.obj | --selftest\n");
    return 2;
}
