/* mesh.h -- attribute-bank mesh for the redundancy experiment.
 *
 * One mesh: nv vertices, nf faces, four channels.  Each channel has its own
 * 16-bit index stream and its own bank of slots.  The position bank has
 * exactly nv slots.  Each attribute bank has more slots than its stream
 * references: the attribute records are deduplicated onto a coarser grid,
 * and an earlier export left the old records in place, so the file carries
 * stale slots that no index ever reaches.
 *
 * Attribute slot ids are a shuffled subset of [0, bank), so they are NOT
 * vertex ids.  Counting a channel through the position stream is a category
 * error; the experiment runs both criteria and compares them.
 */
#ifndef MESH_H
#define MESH_H

#include <stdint.h>
#include <stddef.h>

typedef struct { float x, y, z; } V3;
enum { CH_POS = 0, CH_NRM = 1, CH_UV = 2, CH_COL = 3, CH_N = 4 };
extern const char* const CH_NAME[CH_N];
extern const int CH_COMP[CH_N];

typedef struct {
    const char* name;
    int nv, nf;
    int bank[CH_N];           /* slots stored per channel */
    int live_n[CH_N];         /* slots referenced by the channel's stream */
    V3*   pos;
    V3*   nrm;
    float* uv;
    float* col;
    uint16_t* idx[CH_N];
} Mesh;

/* div_x/div_y: dedup factor of the attribute record grid (1 = full res).
 * bank[ch]: slots actually written to the file (pos must equal nv).
 * Dead slots are filled with drifted copies of live records. */
Mesh  mesh_make(const char* name, int w, int h,
                const int div_x[CH_N], const int div_y[CH_N],
                const int bank[CH_N], float stale_eps, unsigned seed);
void  mesh_free(Mesh* m);

/* counts[s], s in [0,bank).  via_vertex=0: occurrences in the channel's own
 * index stream.  via_vertex=1: the defective rule -- vertex usage of id i
 * applied to attribute row i, for i < min(nv,bank); rows >= nv are never
 * probed and therefore kept. */
void  mesh_counts(const Mesh* m, int ch, int* counts, int via_vertex);

/* keep slots with live[ch][s] != 0, compact banks, remap streams */
Mesh  mesh_prune(const Mesh* src, int* const live[CH_N]);

size_t mesh_bank_bytes(const Mesh* m, int ch);
size_t mesh_attr_bytes(const Mesh* m);
size_t mesh_index_bytes(const Mesh* m);

int    mesh_write(const Mesh* m, const char* path);
Mesh   mesh_read(const char* path);

void   mesh_render(const Mesh* m, int W, int H, float* rgb);
int    pixels_differ(const float* a, const float* b, int W, int H,
                     long* ndiff, float* max_delta);
int    attrs_equal(const Mesh* a, const Mesh* b);
int    write_ppm(const char* path, const float* rgb, int W, int H);
int    write_ppm_diff(const char* path, const float* a, const float* b, int W, int H);

#endif
