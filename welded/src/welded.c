
/* welded: vertex welding by position quantization, with face remap. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 64
#define MAXF 64

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

/* weld vertices closer than eps, remap faces, drop degenerate */
static int mesh_weld(Mesh* m, float eps) {
    if (eps <= 0.0f) return -1;
    V3 uniq[MAXV]; int nu = 0;
    for (int i = 0; i < m->nv; i++) {
        int found = -1;
        for (int j = 0; j < nu; j++) {
            float dx = m->v[i].x - uniq[j].x;
            float dy = m->v[i].y - uniq[j].y;
            float dz = m->v[i].z - uniq[j].z;
            if (dx*dx + dy*dy + dz*dz <= eps*eps) { found = j; break; }
        }
        if (found < 0) { uniq[nu] = m->v[i]; found = nu++; }
        m->map[i] = found;
    }
    int nf = 0;
    for (int i = 0; i < m->nf; i++) {
        int a = m->map[m->f[i].a], b = m->map[m->f[i].b], c = m->map[m->f[i].c];
        if (a == b || b == c || a == c) continue;
        m->f[nf].a = a; m->f[nf].b = b; m->f[nf].c = c; nf++;
    }
    m->nf = nf;
    memcpy(m->v, uniq, sizeof(V3) * nu);
    m->nv = nu;
    m->welded_count = nu;
    return nu;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Mesh m;
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

    check("bad_eps", mesh_weld(&m, 0.0f) == -1);

    Mesh e;
    mesh_init(&e);
    mesh_add_vertex(&e, 0,0,0); mesh_add_vertex(&e, 1,0,0); mesh_add_vertex(&e, 0,1,0);
    mesh_add_tri(&e, 0, 1, 2);
    check("no_weld_needed", mesh_weld(&e, 1e-4f) == 3 && e.nf == 1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
