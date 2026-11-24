#ifndef MML_NORMAL_H
#define MML_NORMAL_H

#include "../../contract/mml_contract.h"

typedef struct mml_normal_options {
    mml_header header;

    /*
     * Accumulator length at or below this value is treated as invalid.
     * Zero means "use the module default". Added as an optional tail field;
     * callers whose header.size stops before it get the default.
     */
    float zero_epsilon;
} mml_normal_options;

#define MML_NORMAL_OPTIONS_V1_SIZE sizeof(mml_header)
#define MML_NORMAL_DEFAULT_ZERO_EPSILON 1e-12f

typedef struct mml_normal_report {
    mml_header header;
    uint32_t zero_accumulator_vertices;
} mml_normal_report;

#define MML_NORMAL_REPORT_V1_SIZE                                              \
    MML_FIELD_END(mml_normal_report, zero_accumulator_vertices)

mml_status mml_normal_options_init(mml_normal_options *options);

/*
 * Compute area-weighted vertex normals using cross(b-a, c-a).
 * out_normals must have out_capacity entries, at least mesh->vertex_count.
 * A vertex whose accumulator is zero/degenerate receives (0,0,0), and the
 * count of such vertices is returned in report.
 */
mml_status mml_compute_vertex_normals(const mml_mesh_view *mesh,
                                      const mml_normal_options *options,
                                      mml_vec3 *out_normals,
                                      uint32_t out_capacity,
                                      mml_normal_report *report);

#endif /* MML_NORMAL_H */
