
/* normal: per-face and area-weighted per-vertex normals. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 32
#define MAXF 32

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

static V3 sub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
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
        out[i].x /= l; out[i].y /= l; out[i].z /= l;
    }
    return bad;   /* number of vertices with zero-length accumulator */
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
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

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
