
/* normal: per-face and area-weighted per-vertex normals.
 *
 * The area-weighted vertex normal computation is split into four units,
 * U1..U4.  Dependencies point one way only:
 *
 *   U1 face_normal          (geometry primitives only)
 *   U2 accumulate_face      (calls U1)
 *   U3 normalize_vertex     (geometry primitives only)
 *   U4 vertex_normals       (orchestrates U2 + U3; the only entry point)
 *
 * No unit calls a unit above it and none call each other sideways, so the
 * graph is a DAG:  U1 <- U2 <- U4,  U3 <- U4.
 *
 * Bit-identity rule: every floating-point expression below is kept in the
 * same order as the original monolith.  In particular accumulation stays a
 * single in-place scatter over `out` in face/vertex order — no temporary
 * per-face buffers are introduced, so the float add sequence is unchanged.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "normal.h"

static V3 sub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 cross(V3 p, V3 q) { V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r; }
static float dot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float len(V3 p) { return sqrtf(dot(p,p)); }
static V3 norm(V3 p) { float l = len(p); V3 r = {0,0,0}; if (l > 0) { r.x=p.x/l; r.y=p.y/l; r.z=p.z/l; } return r; }

static int add_v(Mesh* m, float x, float y, float z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv].x = x; m->v[m->nv].y = y; m->v[m->nv].z = z; return m->nv++;
}
static int add_f(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a<0||b<0||c<0||a>=m->nv||b>=m->nv||c>=m->nv) return -1;
    m->f[m->nf].a=a; m->f[m->nf].b=b; m->f[m->nf].c=c; return m->nf++;
}

/* ---- U1: per-face area-weighted normal --------------------------------
 * Input : the mesh and one face index i.
 * Output: cross(b-a, c-a); its length equals 2*triangle area. Pure, no
 *         mutation. Indices are trusted valid (U4 owns the nf bound). */
V3 face_normal(const Mesh* m, int i) {
    V3 a = m->v[m->f[i].a], b = m->v[m->f[i].b], c = m->v[m->f[i].c];
    return cross(sub(b, a), sub(c, a));
}

/* ---- U2: scatter one face normal onto its three corner accumulators ----
 * Input : mesh, face index i, the accumulator array out[0..nv-1].
 * Output: none by value; out[a], out[b], out[c] are incremented in place.
 *         Same in-place array, same per-component add order as the original
 *         code — no temporary buffer, hence identical rounding. */
void accumulate_face(const Mesh* m, int i, V3 n, V3* out) {
    int ia = m->f[i].a, ib = m->f[i].b, ic = m->f[i].c;
    out[ia].x += n.x; out[ia].y += n.y; out[ia].z += n.z;
    out[ib].x += n.x; out[ib].y += n.y; out[ib].z += n.z;
    out[ic].x += n.x; out[ic].y += n.y; out[ic].z += n.z;
}

/* ---- U3: normalize one accumulated vertex vector -----------------------
 * Input : one accumulator V3 (in/out).
 * Output: 0 if len < 1e-12, in which case the vector is left untouched
 *         (matching the original `continue`); otherwise 1 and *v is divided
 *         componentwise by len. The 1e-12 threshold is deliberately not
 *         `> 0`: note it differs from the unused norm() helper above. */
int normalize_vertex(V3* v) {
    float l = len(*v);
    if (l < 1e-12f) return 0;
    v->x /= l; v->y /= l; v->z /= l;
    return 1;
}

/* ---- U4: orchestrator --------------------------------------------------
 * Input : mesh m, caller-provided out[0..nv-1].
 * Output: -1 if the mesh has no faces; otherwise the count of vertices whose
 *         accumulator was below the zero-length threshold. Owns the empty
 *         guard, zero initialization, and both loop orders. */
int vertex_normals(const Mesh* m, V3* out) {
    if (m->nf == 0) return -1;
    for (int i = 0; i < m->nv; i++) { out[i].x = out[i].y = out[i].z = 0.0f; }
    for (int i = 0; i < m->nf; i++) {
        V3 n = face_normal(m, i);       /* U1 */
        accumulate_face(m, i, n, out);  /* U2 */
    }
    int bad = 0;
    for (int i = 0; i < m->nv; i++) {
        if (!normalize_vertex(&out[i])) bad++;  /* U3 */
    }
    return bad;
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
