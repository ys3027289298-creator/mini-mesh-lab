/* welded: contract-based vertex welding, face remap and degenerate removal. */
#include "../include/mml_welded.h"

#include <stdio.h>
#include <stdlib.h>

mml_status mml_weld_options_init(mml_weld_options *options, float epsilon)
{
    if (options == NULL)
        return MML_E_NULL;
    MML_STRUCT_INIT(options, mml_weld_options);
    options->epsilon = epsilon;
    return MML_OK;
}

mml_status mml_weld_mesh(const mml_mesh_view *in,
                         const mml_weld_options *options,
                         mml_mesh_buffer *out,
                         mml_weld_report *report)
{
    mml_status st;
    mml_index *owned_map = NULL;
    mml_index *map;
    uint32_t i;
    uint32_t dropped = 0;

    if (in == NULL || options == NULL || out == NULL)
        return MML_E_NULL;

    st = mml_mesh_view_validate(in);
    if (st != MML_OK)
        return st;

    st = mml_check_header(&options->header, MML_WELD_OPTIONS_V1_SIZE);
    if (st != MML_OK)
        return st;
    if (!isfinite(options->epsilon) || options->epsilon <= 0.0f)
        return MML_E_ARGUMENT;

    st = mml_check_header(&out->header, MML_MESH_BUFFER_V1_SIZE);
    if (st != MML_OK)
        return st;

    if (report != NULL) {
        st = mml_check_header(&report->header, MML_WELD_REPORT_V1_SIZE);
        if (st != MML_OK)
            return st;
    }

    /* Welding never creates more vertices/triangles than its input. */
    if (in->vertex_count > 0) {
        if (out->positions == NULL ||
            out->vertex_capacity < in->vertex_count)
            return MML_E_CAPACITY;
    }
    if (in->triangle_count > 0) {
        if (out->triangles == NULL ||
            out->triangle_capacity < in->triangle_count)
            return MML_E_CAPACITY;
    }

    map = NULL;
    if (options->old_to_new != NULL) {
        if (options->old_to_new_capacity < in->vertex_count)
            return MML_E_CAPACITY;
        map = options->old_to_new;
    }

    if (map == NULL && in->vertex_count > 0) {
        map = calloc(in->vertex_count, sizeof(*map));
        if (map == NULL)
            return MML_E_NOMEM;
        owned_map = map;
    }

    out->vertex_count = 0;
    out->triangle_count = 0;
    /* Welding invalidates derived channels such as normals. */
    out->flags = 0;

    for (i = 0; i < in->vertex_count; ++i) {
        mml_index found = MML_INVALID_INDEX;
        uint32_t j;

        for (j = 0; j < out->vertex_count; ++j) {
            double dx = (double)in->positions[i].x - out->positions[j].x;
            double dy = (double)in->positions[i].y - out->positions[j].y;
            double dz = (double)in->positions[i].z - out->positions[j].z;
            double eps = options->epsilon;

            if (dx * dx + dy * dy + dz * dz <= eps * eps) {
                found = j;
                break;
            }
        }

        if (found == MML_INVALID_INDEX) {
            /* Capacity was checked against the worst case above. */
            found = out->vertex_count++;
            out->positions[found] = in->positions[i];
        }
        if (map != NULL)
            map[i] = found;
    }

    for (i = 0; i < in->triangle_count; ++i) {
        const mml_triangle *src = &in->triangles[i];
        mml_index a = map[src->a];
        mml_index b = map[src->b];
        mml_index c = map[src->c];

        if (a == b || b == c || a == c) {
            ++dropped;
            continue;
        }

        /* Capacity was checked against the worst case above. */
        out->triangles[out->triangle_count++] =
            (mml_triangle){ a, b, c };
    }

    if (report != NULL) {
        report->vertices_before = in->vertex_count;
        report->vertices_after = out->vertex_count;
        report->triangles_before = in->triangle_count;
        report->triangles_after = out->triangle_count;
        report->degenerate_triangles_dropped = dropped;
    }

    free(owned_map);
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

static mml_status add_tri(mml_triangle *tris,
                          uint32_t capacity,
                          uint32_t *count,
                          uint32_t vertex_count,
                          mml_index a,
                          mml_index b,
                          mml_index c)
{
    if (*count >= capacity)
        return MML_E_CAPACITY;
    if (a >= vertex_count || b >= vertex_count || c >= vertex_count)
        return MML_E_INDEX;
    if (a == b || b == c || a == c)
        return MML_E_ARGUMENT;
    tris[(*count)++] = (mml_triangle){ a, b, c };
    return MML_OK;
}

int main(void)
{
    mml_vec3 positions[5] = {
        { 0.0f, 0.0f, 0.0f },
        { 1.0f, 0.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { 0.0f, 0.0f, 0.0f },        /* duplicate of 0 */
        { 1.0001f, 0.0f, 0.0f }      /* near duplicate of 1 */
    };
    mml_triangle tris[4];
    uint32_t triangle_count = 0;

    check("added",
          sizeof(positions) / sizeof(positions[0]) == 5);
    check("tri_ok",
          add_tri(tris, 4, &triangle_count, 5, 0, 1, 2) == MML_OK);
    check("tri_dup_rejected",
          add_tri(tris, 4, &triangle_count, 5, 0, 0, 2) < 0);
    check("tri_range_rejected",
          add_tri(tris, 4, &triangle_count, 5, 0, 1, 99) == MML_E_INDEX);
    add_tri(tris, 4, &triangle_count, 5, 3, 4, 2);
    add_tri(tris, 4, &triangle_count, 5, 0, 3, 1);

    mml_mesh_view in;
    MML_STRUCT_INIT(&in, mml_mesh_view);
    in.vertex_count = 5;
    in.triangle_count = triangle_count;
    in.positions = positions;
    in.triangles = tris;

    mml_vec3 welded_positions[5];
    mml_triangle welded_tris[4];
    mml_mesh_buffer out;
    MML_STRUCT_INIT(&out, mml_mesh_buffer);
    out.vertex_capacity = 5;
    out.triangle_capacity = 4;
    out.positions = welded_positions;
    out.triangles = welded_tris;

    mml_index map[5];
    mml_weld_options options;
    mml_weld_options_init(&options, 1e-3f);
    options.old_to_new = map;
    options.old_to_new_capacity = 5;

    mml_weld_report report;
    MML_STRUCT_INIT(&report, mml_weld_report);

    check("weld_ok", mml_weld_mesh(&in, &options, &out, &report) == MML_OK);
    check("weld_count", report.vertices_after == 3);
    check("map_zero", map[0] == map[3]);
    check("map_one", map[1] == map[4]);
    check("degenerate_dropped", report.triangles_after == 2);

    mml_weld_options bad;
    mml_weld_options_init(&bad, 0.0f);
    check("bad_eps", mml_weld_mesh(&in, &bad, &out, NULL) < 0);

    mml_vec3 p[3] = {
        { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }
    };
    mml_triangle f[1] = { { 0, 1, 2 } };
    mml_mesh_view e;
    MML_STRUCT_INIT(&e, mml_mesh_view);
    e.vertex_count = 3;
    e.triangle_count = 1;
    e.positions = p;
    e.triangles = f;

    mml_vec3 ep[3];
    mml_triangle ef[1];
    mml_mesh_buffer eo;
    MML_STRUCT_INIT(&eo, mml_mesh_buffer);
    eo.vertex_capacity = 3;
    eo.triangle_capacity = 1;
    eo.positions = ep;
    eo.triangles = ef;

    mml_weld_report er;
    MML_STRUCT_INIT(&er, mml_weld_report);
    mml_weld_options eo_opt;
    mml_weld_options_init(&eo_opt, 1e-4f);

    check("no_weld_needed",
          mml_weld_mesh(&e, &eo_opt, &eo, &er) == MML_OK &&
          er.vertices_after == 3 && er.triangles_after == 1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
