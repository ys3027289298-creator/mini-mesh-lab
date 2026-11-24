/* bbox: contract-based AABB and sliver detection. */
#include "../include/mml_bbox.h"

#include <stdio.h>

mml_status mml_bbox_options_init(mml_bbox_options *options)
{
    if (options == NULL)
        return MML_E_NULL;
    MML_STRUCT_INIT(options, mml_bbox_options);
    return MML_OK;
}

mml_status mml_bbox_compute(const mml_mesh_view *mesh,
                            const mml_bbox_options *options,
                            mml_box *box)
{
    mml_status st;
    uint32_t i;

    if (mesh == NULL || box == NULL)
        return MML_E_NULL;

    st = mml_mesh_view_validate(mesh);
    if (st != MML_OK)
        return st;

    if (options != NULL) {
        st = mml_check_header(&options->header,
                              MML_BBOX_OPTIONS_V1_SIZE);
        if (st != MML_OK)
            return st;
    }

    st = mml_check_header(&box->header, MML_BOX_V1_SIZE);
    if (st != MML_OK)
        return st;

    if (mesh->vertex_count == 0) {
        box->min = (mml_vec3){ 0.0f, 0.0f, 0.0f };
        box->max = (mml_vec3){ 0.0f, 0.0f, 0.0f };
        box->valid = 0;
        return MML_E_EMPTY;
    }

    box->min = mesh->positions[0];
    box->max = mesh->positions[0];
    for (i = 1; i < mesh->vertex_count; ++i) {
        mml_vec3 p = mesh->positions[i];

        if (p.x < box->min.x) box->min.x = p.x;
        if (p.y < box->min.y) box->min.y = p.y;
        if (p.z < box->min.z) box->min.z = p.z;
        if (p.x > box->max.x) box->max.x = p.x;
        if (p.y > box->max.y) box->max.y = p.y;
        if (p.z > box->max.z) box->max.z = p.z;
    }
    box->valid = 1;
    return MML_OK;
}

mml_status mml_box_size(const mml_box *box, mml_vec3 *size)
{
    mml_status st;

    if (box == NULL || size == NULL)
        return MML_E_NULL;
    st = mml_check_header(&box->header, MML_BOX_V1_SIZE);
    if (st != MML_OK)
        return st;
    if (!box->valid)
        return MML_E_EMPTY;

    size->x = box->max.x - box->min.x;
    size->y = box->max.y - box->min.y;
    size->z = box->max.z - box->min.z;
    return MML_OK;
}

mml_status mml_box_diagonal(const mml_box *box, float *diagonal)
{
    mml_vec3 size;
    mml_status st;

    if (diagonal == NULL)
        return MML_E_NULL;
    st = mml_box_size(box, &size);
    if (st != MML_OK)
        return st;

    *diagonal = sqrtf(size.x * size.x + size.y * size.y +
                      size.z * size.z);
    return MML_OK;
}

mml_status mml_box_is_sliver(const mml_box *box,
                             float ratio,
                             int *is_sliver)
{
    mml_vec3 size;
    float minimum_extent;
    float diagonal;
    mml_status st;

    if (box == NULL || is_sliver == NULL)
        return MML_E_NULL;
    st = mml_check_header(&box->header, MML_BOX_V1_SIZE);
    if (st != MML_OK)
        return st;
    if (!box->valid)
        return MML_E_EMPTY;
    if (!isfinite(ratio) || ratio <= 0.0f || ratio >= 1.0f)
        return MML_E_ARGUMENT;

    st = mml_box_size(box, &size);
    if (st != MML_OK)
        return st;

    minimum_extent = size.x < size.y ? size.x : size.y;
    if (size.z < minimum_extent)
        minimum_extent = size.z;

    st = mml_box_diagonal(box, &diagonal);
    if (st != MML_OK)
        return st;
    if (diagonal <= 0.0f)
        return MML_E_DEGENERATE;

    *is_sliver = (minimum_extent < ratio * diagonal);
    return MML_OK;
}

#ifndef MML_NO_SELFTEST
static int failures = 0;

static void check(const char *name, int ok)
{
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok)
        ++failures;
}

int main(void)
{
    mml_box b;
    mml_box_init(&b);
    check("empty_invalid", b.valid == 0);
    int sliver = -1;
    check("empty_sliver_err",
          mml_box_is_sliver(&b, 0.1f, &sliver) == MML_E_EMPTY);

    mml_vec3 one_p[1] = { { 0, 0, 0 } };
    mml_mesh_view one;
    MML_STRUCT_INIT(&one, mml_mesh_view);
    one.vertex_count = 1;
    one.positions = one_p;
    check("single_valid", mml_bbox_compute(&one, NULL, &b) == MML_OK &&
                          b.valid == 1);
    float diagonal = -1.0f;
    check("single_zero_diag",
          mml_box_diagonal(&b, &diagonal) == MML_OK &&
          fabsf(diagonal) < 1e-6f);

    mml_vec3 cube_p[2] = { { 0, 0, 0 }, { 1, 1, 1 } };
    mml_mesh_view cube;
    MML_STRUCT_INIT(&cube, mml_mesh_view);
    cube.vertex_count = 2;
    cube.positions = cube_p;
    mml_box_init(&b);
    mml_bbox_compute(&cube, NULL, &b);
    mml_vec3 s;
    mml_box_size(&b, &s);
    check("cube_size",
          fabsf(s.x - 1.0f) < 1e-6f &&
          fabsf(s.y - 1.0f) < 1e-6f &&
          fabsf(s.z - 1.0f) < 1e-6f);
    check("cube_not_sliver",
          mml_box_is_sliver(&b, 0.1f, &sliver) == MML_OK &&
          sliver == 0);

    mml_vec3 sl_p[2] = { { 0, 0, 0 }, { 10, 10, 0.01f } };
    mml_mesh_view sl_mesh;
    MML_STRUCT_INIT(&sl_mesh, mml_mesh_view);
    sl_mesh.vertex_count = 2;
    sl_mesh.positions = sl_p;
    mml_box sl;
    mml_box_init(&sl);
    mml_bbox_compute(&sl_mesh, NULL, &sl);
    check("sliver_detected",
          mml_box_is_sliver(&sl, 0.1f, &sliver) == MML_OK &&
          sliver == 1);

    check("bad_ratio",
          mml_box_is_sliver(&sl, 0.0f, &sliver) == MML_E_ARGUMENT &&
          mml_box_is_sliver(&sl, 1.0f, &sliver) == MML_E_ARGUMENT);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
