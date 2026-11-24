/*
 * Minimal end-to-end contract scenario:
 *
 *   raw local mesh -> weld -> contract mesh view -> normals
 *                  -> same contract view (+normals) -> bbox/sliver
 *
 * There is no vertex-layout conversion between modules: every boundary uses
 * mml_vec3 / mml_triangle / mml_mesh_view.
 */
#include "mml_contract.h"
#include "mml_welded.h"
#include "mml_normal.h"
#include "mml_bbox.h"

#include <math.h>
#include <stdio.h>

static int failures = 0;

static void expect(const char *name, int ok)
{
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok)
        ++failures;
}

static int near_float(float a, float b)
{
    return fabsf(a - b) < 1e-5f;
}

int main(void)
{
    /* Two exact/near duplicates plus one triangle that becomes degenerate. */
    mml_vec3 input_positions[5] = {
        { 0.0f, 0.0f, 0.0f },
        { 1.0f, 0.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { 0.0f, 0.0f, 0.0f },
        { 1.0001f, 0.0f, 0.0f }
    };
    mml_triangle input_triangles[3] = {
        { 0, 1, 2 }, /* canonical triangle */
        { 3, 4, 2 }, /* duplicate after welding */
        { 0, 3, 1 }  /* collapses to (0,0,1) and is dropped */
    };

    mml_mesh_view input;
    MML_STRUCT_INIT(&input, mml_mesh_view);
    input.vertex_count = 5;
    input.triangle_count = 3;
    input.positions = input_positions;
    input.triangles = input_triangles;

    /* Caller-owned weld output storage. */
    mml_vec3 welded_positions[5];
    mml_triangle welded_triangles[3];
    mml_mesh_buffer welded_buffer;
    MML_STRUCT_INIT(&welded_buffer, mml_mesh_buffer);
    welded_buffer.vertex_capacity = 5;
    welded_buffer.triangle_capacity = 3;
    welded_buffer.positions = welded_positions;
    welded_buffer.triangles = welded_triangles;

    mml_index old_to_new[5];
    mml_weld_options weld_options;
    mml_weld_options_init(&weld_options, 1e-3f);
    weld_options.old_to_new = old_to_new;
    weld_options.old_to_new_capacity = 5;

    mml_weld_report weld_report;
    MML_STRUCT_INIT(&weld_report, mml_weld_report);

    mml_status st = mml_weld_mesh(&input, &weld_options,
                                  &welded_buffer, &weld_report);
    expect("weld_status_ok", st == MML_OK);
    expect("welded_three_vertices", weld_report.vertices_after == 3);
    expect("welded_two_triangles", weld_report.triangles_after == 2);
    expect("weld_dropped_one",
           weld_report.degenerate_triangles_dropped == 1);
    expect("weld_local_map",
           old_to_new[0] == 0 && old_to_new[3] == 0 &&
           old_to_new[1] == 1 && old_to_new[4] == 1 &&
           old_to_new[2] == 2);

    printf("weld: %u -> %u vertices, %u -> %u triangles\n",
           weld_report.vertices_before, weld_report.vertices_after,
           weld_report.triangles_before, weld_report.triangles_after);
    printf("map(local old -> local new):");
    for (uint32_t i = 0; i < 5; ++i)
        printf(" %u->%u", i, old_to_new[i]);
    printf("\n");

    /* Borrow the weld output as a read-only contract view. */
    mml_mesh_view welded_view;
    MML_STRUCT_INIT(&welded_view, mml_mesh_view);
    welded_view.vertex_count = welded_buffer.vertex_count;
    welded_view.triangle_count = welded_buffer.triangle_count;
    welded_view.positions = welded_buffer.positions;
    welded_view.triangles = welded_buffer.triangles;

    mml_vec3 normals[5];
    mml_normal_report normal_report;
    MML_STRUCT_INIT(&normal_report, mml_normal_report);
    st = mml_compute_vertex_normals(&welded_view, NULL, normals, 5,
                                    &normal_report);
    expect("normal_status_ok", st == MML_OK);
    expect("no_invalid_normals",
           normal_report.zero_accumulator_vertices == 0);
    expect("normal_points_positive_z",
           near_float(normals[0].x, 0.0f) &&
           near_float(normals[0].y, 0.0f) &&
           near_float(normals[0].z, 1.0f));

    /* The same view can now advertise the caller-owned derived channel. */
    welded_view.vertex_normals = normals;
    welded_view.flags |= MML_MESH_HAS_NORMALS;

    printf("normal[0]=(%.6f %.6f %.6f)\n",
           normals[0].x, normals[0].y, normals[0].z);

    mml_box box;
    mml_box_init(&box);
    st = mml_bbox_compute(&welded_view, NULL, &box);
    expect("bbox_status_ok", st == MML_OK);
    expect("bbox_min",
           near_float(box.min.x, 0.0f) &&
           near_float(box.min.y, 0.0f) &&
           near_float(box.min.z, 0.0f));
    expect("bbox_max",
           near_float(box.max.x, 1.0f) &&
           near_float(box.max.y, 1.0f) &&
           near_float(box.max.z, 0.0f));

    int sliver = 0;
    st = mml_box_is_sliver(&box, 0.1f, &sliver);
    expect("bbox_sliver_status_ok", st == MML_OK);
    expect("flat_triangle_is_sliver", sliver == 1);

    printf("bbox min=(%.3f %.3f %.3f) max=(%.3f %.3f %.3f) sliver=%d\n",
           box.min.x, box.min.y, box.min.z,
           box.max.x, box.max.y, box.max.z, sliver);

    printf("chain_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
