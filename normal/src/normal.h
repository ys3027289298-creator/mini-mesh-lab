#ifndef NORMAL_H
#define NORMAL_H

#define MAXV 32
#define MAXF 32

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { V3 v[MAXV]; int nv; Tri f[MAXF]; int nf; } Mesh;

/* U4: the only externally visible entry point of the split pipeline. */
int vertex_normals(const Mesh* m, V3* out);

#endif
