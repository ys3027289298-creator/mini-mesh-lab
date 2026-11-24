/* seam: tolerance-driven normal merging across UV seams.
 *
 * Split vertex records (same world position, separate UV island) get independent
 * per-record normals and show a shading line at every UV seam.  This module welds
 * records by position tolerance epsilon, optionally with a normal-cone guard that
 * keeps genuine crease edges from merging, then shares one normal per welded
 * group (the "open UV channel" step). */
#ifndef SEAM_H
#define SEAM_H

#define MAXV 4096
#define MAXF 8192

/* folded-strip generator parameters, exposed so tests can predict the geometry */
#define SEAM_PANELS        4
#define SEAM_COL_PER_PANEL 8
#define SEAM_ROW           24
#define SEAM_CELL          (1.0f / SEAM_COL_PER_PANEL)
#define SEAM_F_SEAM        0.02f   /* seam twin gap / cell */
#define SEAM_F_CREASE      0.20f   /* crease neighbour gap / cell */

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct { float u, v; } UV;

typedef struct {
    V3  pos[MAXV];      /* split vertex records */
    V3  nrm[MAXV];      /* per-record area-weighted normal (filled by mesh_normals) */
    V3  truth[MAXV];    /* analytic smooth normal at the record's position */
    UV  uv[MAXV];
    int seam[MAXV];     /* 1 if this record was duplicated by a UV seam split */
    int nv;
    Tri f[MAXF];
    int nf;
} Mesh;

typedef struct {
    int a, b;           /* two record indices */
    float gap;          /* world-space distance between them */
} Pair;

#define MAXPAIRS 4096

typedef struct {
    Pair seam[MAXPAIRS];   int nseam;    /* same-panel UV-seam twins (should merge) */
    Pair crease[MAXPAIRS]; int ncrease;  /* fold-line neighbours (must never merge) */
} PairSets;

void mesh_init(Mesh* m);
int  mesh_add_vertex(Mesh* m, V3 p, UV uv, int seam_flag);
int  mesh_add_tri(Mesh* m, int a, int b, int c);

/* folded-strip test geometry: P flat tilted panels separated by physical crease
 * gaps, with nseam internal UV splits distributed across the panels.
 * bump amplitude and all distances scale with `scale`. */
void mesh_build_folded_strip(Mesh* m, int nseam, float scale);

/* median edge length of the triangulation: the patch-size measure L */
float mesh_characteristic_length(const Mesh* m);
void  mesh_normals(Mesh* m);
void  mesh_ground_truth_normals(const Mesh* m, V3* out);

/* classify split records into seam twins and crease neighbours */
void  mesh_find_pairs(const Mesh* m, PairSets* ps);

/* union-find over vertex records */
typedef struct { int parent[MAXV]; int n; } UF;
void  uf_init(UF* uf, int n);
int   uf_find(UF* uf, int x);
void  uf_union(UF* uf, int a, int b);

/* position-only weld: every record pair closer than eps is joined */
void  weld_position(const Mesh* m, float eps, UF* uf);
/* guarded weld: same, but a pair whose per-record normals differ by more than
 * angle_limit_deg is refused (crease protection) */
void  weld_guarded(const Mesh* m, float eps, float angle_limit_deg, UF* uf);

/* one averaged normal per weld group, written back to every member record */
void  share_group_normals(const Mesh* m, const UF* uf, V3* out);

/* number of connected UV islands (records linked by shared triangle edges) */
int   count_uv_channels(const Mesh* m);
/* number of normal groups that still split a seam position after welding */
int   count_open_normal_channels(const Mesh* m, const UF* uf, const PairSets* ps);

/* metrics over a finished normal field `out` */
int   count_unmerged_seams(const PairSets* ps, const UF* uf);
int   count_mismerged_creases(const PairSets* ps, const UF* uf);
/* seam pairs whose two records still carry different normals (bit-exact):
 * 0 once a normal is shared through the weld group, every pair otherwise */
int   count_open_normal_seams(const PairSets* ps, const V3* nrm);
/* distinct normal values held inside each weld group, summed over groups:
 * one per group after channel opening, one per record before */
int   count_normal_values(const Mesh* m, const UF* uf, const V3* nrm);
float mean_normal_error_deg(const Mesh* m, const V3* out);

/* orthographic software rasteriser.
 * renders shading from `nrm` and from m->truth; returns pixels inside seam zones
 * whose shading differs from the truth render by more than shade_tol.
 * when ppm_path is non-NULL the test render is written there (binary PPM). */
#define RENDER_W 960
#define RENDER_H 520
int   render_bright_pixels(const Mesh* m, const V3* nrm, float shade_tol,
                           const char* ppm_path);

#endif
