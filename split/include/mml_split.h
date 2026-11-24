/* split: cut a triangle mesh into connected components.
 *
 * Each component gets its own compact vertex array and its faces are
 * renumbered against that array, so every component's indices start at 0.
 * Vertices are never shared between components, and the totals are
 * conserved: sum(nv) == input nv, sum(nf) == input nf (no face is dropped,
 * because splitting happens before any welding).
 */
#ifndef MML_SPLIT_H
#define MML_SPLIT_H

#include <stddef.h>

#define MML_SPLIT_MAXV 256
#define MML_SPLIT_MAXF 256

typedef struct { float x, y, z; } MmlV3;
typedef struct { int a, b, c; } MmlTri;

typedef struct {
    MmlV3 v[MML_SPLIT_MAXV]; int nv;
    MmlTri f[MML_SPLIT_MAXF]; int nf;
} MmlMesh;

typedef struct {
    MmlMesh m;
    int src[MML_SPLIT_MAXV];   /* local index -> source vertex index */
} MmlComp;

/* Split into connected components. Returns the component count, or -1 if
 * the mesh exceeds the fixed capacities or out[] is too small.
 * Components are ordered by the smallest source vertex index they own.
 * Isolated vertices (not referenced by any face) become one-vertex
 * components with no faces. */
int mml_split(const MmlMesh* in, MmlComp* out, int maxc);

/* Integrity checks, each returns 1 when the property holds. */
int mml_check_in_range(const MmlComp* c, int nc);
int mml_check_dense(const MmlComp* c, int nc);       /* every local vertex used */
int mml_check_disjoint(const MmlComp* c, int nc);    /* no shared source vertex */
int mml_check_conserved(const MmlMesh* in, const MmlComp* c, int nc);
int mml_check_faces(const MmlMesh* in, const MmlComp* c, int nc);

#endif
