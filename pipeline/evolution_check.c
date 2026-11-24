/*
 * Backward/forward size-header evolution check.
 *
 * It deliberately presents a shorter "old caller" options header and a
 * larger "future caller" size to verify that callees use size gating rather
 * than exact struct-size equality.
 */
#include "mml_contract.h"
#include "mml_normal.h"

#include <math.h>
#include <stdio.h>

static int failures = 0;

static void expect(const char *name, int ok)
{
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok)
        ++failures;
}

int main(void)
{
    mml_vec3 positions[3] = {
        { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }
    };
    mml_triangle triangles[1] = { { 0, 1, 2 } };
    mml_mesh_view mesh;
    MML_STRUCT_INIT(&mesh, mml_mesh_view);
    mesh.vertex_count = 3;
    mesh.triangle_count = 1;
    mesh.positions = positions;
    mesh.triangles = triangles;

    mml_vec3 normals[3];
    mml_normal_report report;

    /* Old caller: header says it only knows the mml_header prefix. */
    mml_normal_options old_options;
    mml_normal_options_init(&old_options);
    old_options.zero_epsilon = 1.5f; /* must be ignored due to short size */
    old_options.header.size = MML_NORMAL_OPTIONS_V1_SIZE;

    MML_STRUCT_INIT(&report, mml_normal_report);
    mml_status st = mml_compute_vertex_normals(
        &mesh, &old_options, normals, 3, &report);
    expect("short_old_size_accepted", st == MML_OK);
    expect("short_old_size_uses_default",
           report.zero_accumulator_vertices == 0 &&
           fabsf(fabsf(normals[0].z) - 1.0f) < 1e-5f);

    /* Current caller exposes the tail field; 1.5 marks length-1 accumulators. */
    mml_normal_options current_options;
    mml_normal_options_init(&current_options);
    current_options.zero_epsilon = 1.5f;

    MML_STRUCT_INIT(&report, mml_normal_report);
    st = mml_compute_vertex_normals(
        &mesh, &current_options, normals, 3, &report);
    expect("appended_field_is_read",
           st == MML_OK &&
           report.zero_accumulator_vertices == 3 &&
           normals[0].x == 0.0f && normals[0].y == 0.0f &&
           normals[0].z == 0.0f);

    /* A newer caller may have appended more bytes; callee must not reject. */
    current_options.zero_epsilon = 0.0f; /* zero means default */
    current_options.header.size = 64;    /* deliberately larger than object */

    MML_STRUCT_INIT(&report, mml_normal_report);
    st = mml_compute_vertex_normals(
        &mesh, &current_options, normals, 3, &report);
    expect("larger_future_size_accepted",
           st == MML_OK &&
           report.zero_accumulator_vertices == 0 &&
           fabsf(fabsf(normals[0].z) - 1.0f) < 1e-5f);

    printf("evolution_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
