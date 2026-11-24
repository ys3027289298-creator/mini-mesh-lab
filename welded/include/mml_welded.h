#ifndef MML_WELDED_H
#define MML_WELDED_H

#include "../../contract/mml_contract.h"

/*
 * Welding options. old_to_new, when supplied by a caller built with the full
 * structure, is a caller-owned output array of in->vertex_count entries.
 * MML_INVALID_INDEX means "no output vertex"; the current welder always maps
 * every input vertex to some output vertex.
 */
typedef struct mml_weld_options {
    mml_header header;
    float epsilon;

    /* Optional caller-owned output; NULL lets the module use temporary memory. */
    mml_index *old_to_new;
    uint32_t old_to_new_capacity;
} mml_weld_options;

#define MML_WELD_OPTIONS_V1_SIZE                                              \
    MML_FIELD_END(mml_weld_options, old_to_new_capacity)

typedef struct mml_weld_report {
    mml_header header;
    uint32_t vertices_before;
    uint32_t vertices_after;
    uint32_t triangles_before;
    uint32_t triangles_after;
    uint32_t degenerate_triangles_dropped;
} mml_weld_report;

#define MML_WELD_REPORT_V1_SIZE                                                \
    MML_FIELD_END(mml_weld_report, degenerate_triangles_dropped)

mml_status mml_weld_options_init(mml_weld_options *options, float epsilon);

/*
 * Weld by squared distance <= epsilon^2, remap triangles and drop triangles
 * whose remapped indices repeat. Input and output storage must not overlap.
 */
mml_status mml_weld_mesh(const mml_mesh_view *in,
                         const mml_weld_options *options,
                         mml_mesh_buffer *out,
                         mml_weld_report *report);

#endif /* MML_WELDED_H */
