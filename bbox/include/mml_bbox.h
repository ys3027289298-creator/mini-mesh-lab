#ifndef MML_BBOX_H
#define MML_BBOX_H

#include "../../contract/mml_contract.h"

/* Present so future tuning fields can be appended without changing calls. */
typedef struct mml_bbox_options {
    mml_header header;
} mml_bbox_options;

#define MML_BBOX_OPTIONS_V1_SIZE sizeof(mml_header)

mml_status mml_bbox_options_init(mml_bbox_options *options);

/* Build an AABB. An empty view returns MML_E_EMPTY and sets box->valid=0. */
mml_status mml_bbox_compute(const mml_mesh_view *mesh,
                            const mml_bbox_options *options,
                            mml_box *box);

mml_status mml_box_size(const mml_box *box, mml_vec3 *size);
mml_status mml_box_diagonal(const mml_box *box, float *diagonal);

/*
 * Sliver test: smallest extent < ratio * diagonal.
 * Returns MML_OK and sets *is_sliver to 0/1. An invalid box, invalid ratio
 * or zero diagonal returns a negative mml_status.
 */
mml_status mml_box_is_sliver(const mml_box *box,
                             float ratio,
                             int *is_sliver);

#endif /* MML_BBOX_H */
