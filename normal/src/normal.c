/* normal: area-weighted vertex normals on the shared mesh contract. */
#include "../include/mml_normal.h"

#include <stdio.h>

static mml_vec3 vec_sub(mml_vec3 p, mml_vec3 q)
{
    return (mml_vec3){ p.x - q.x, p.y - q.y, p.z - q.z };
}

static mml_vec3 vec_cross(mml_vec3 p, mml_vec3 q)
{
    return (mml_vec3){
        p.y * q.z - p.z * q.y,
        p.z * q.x - p.x * q.z,
        p.x * q.y - p.y * q.x
    };
}

static float vec_dot(mml_vec3 p, mml_vec3 q)
{
    return p.x * q.x + p.y * q.y + p.z * q.z;
}

static float vec_length(mml_vec3 p)
{
    return sqrtf(vec_dot(p, p));
}

mml_status mml_normal_options_init(mml_normal_options *options)
{
    if (options == NULL)
        return MML_E_NULL;
    MML_STRUCT_INIT(options, mml_normal_options);
    options->zero_epsilon = 0.0f; /* zero selects the documented default */
    return MML_OK;
}

mml_status mml_compute_vertex_normals(const mml_mesh_view *mesh,
                                      const mml_normal_options *options,
                                      mml_vec3 *out_normals,
                                      uint32_t out_capacity,
                                      mml_normal_report *report)
{
    mml_status st;
    float zero_epsilon = MML_NORMAL_DEFAULT_ZERO_EPSILON;
    uint32_t bad = 0;
    uint32_t i;

    if (mesh == NULL || out_normals == NULL)
        return MML_E_NULL;

    st = mml_mesh_view_validate(mesh);
    if (st != MML_OK)
        return st;

    if (options != NULL) {
        st = mml_check_header(&options->header,
                              MML_NORMAL_OPTIONS_V1_SIZE);
        if (st != MML_OK)
            return st;

        /* Size-gated tail field: old callers do not expose zero_epsilon. */
        if (MML_HAS_FIELD(&options->header, mml_normal_options,
                          zero_epsilon) &&
            options->zero_epsilon > 0.0f)
            zero_epsilon = options->zero_epsilon;
    }

    if (report != NULL) {
        st = mml_check_header(&report->header, MML_NORMAL_REPORT_V1_SIZE);
        if (st != MML_OK)
            return st;
    }

    if (out_capacity < mesh->vertex_count)
        return MML_E_CAPACITY;
    if (mesh->triangle_count == 0)
        return MML_E_EMPTY;

    for (i = 0; i < mesh->vertex_count; ++i)
        out_normals[i] = (mml_vec3){ 0.0f, 0.0f, 0.0f };

    /* Accumulate un-normalized cross products: their length is 2*area. */
    for (i = 0; i < mesh->triangle_count; ++i) {
        const mml_triangle *f = &mesh->triangles[i];
        mml_vec3 a = mesh->positions[f->a];
        mml_vec3 b = mesh->positions[f->b];
        mml_vec3 c = mesh->positions[f->c];
        mml_vec3 n = vec_cross(vec_sub(b, a), vec_sub(c, a));

        out_normals[f->a].x += n.x;
        out_normals[f->a].y += n.y;
        out_normals[f->a].z += n.z;
        out_normals[f->b].x += n.x;
        out_normals[f->b].y += n.y;
        out_normals[f->b].z += n.z;
        out_normals[f->c].x += n.x;
        out_normals[f->c].y += n.y;
        out_normals[f->c].z += n.z;
    }

    for (i = 0; i < mesh->vertex_count; ++i) {
        float length = vec_length(out_normals[i]);

        if (!isfinite(length))
            return MML_E_OVERFLOW;

        /* Contract invalid-normal value is the all-zero vector. */
        if (length <= zero_epsilon) {
            out_normals[i] = (mml_vec3){ 0.0f, 0.0f, 0.0f };
            ++bad;
            continue;
        }

        out_normals[i].x /= length;
        out_normals[i].y /= length;
        out_normals[i].z /= length;
    }

    if (report != NULL)
        report->zero_accumulator_vertices = bad;
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
    mml_vec3 p[3] = {
        { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }
    };
    mml_triangle f[1] = { { 0, 1, 2 } };
    mml_mesh_view m;
    MML_STRUCT_INIT(&m, mml_mesh_view);
    m.vertex_count = 3;
    m.triangle_count = 1;
    m.positions = p;
    m.triangles = f;

    mml_vec3 n[3];
    mml_normal_report report;
    MML_STRUCT_INIT(&report, mml_normal_report);
    check("no_degenerate",
          mml_compute_vertex_normals(&m, NULL, n, 3, &report) == MML_OK &&
          report.zero_accumulator_vertices == 0);
    check("unit_length", fabsf(vec_length(n[0]) - 1.0f) < 1e-5f);
    check("points_z", fabsf(fabsf(n[0].z) - 1.0f) < 1e-5f);

    /* Different-area triangles sharing an edge: weighting favours big one. */
    mml_vec3 wp[4] = {
        { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 2 }
    };
    mml_triangle wf[2] = {
        { 0, 1, 2 }, /* small, normal +z */
        { 0, 1, 3 }  /* large, normal -y-ish */
    };
    mml_mesh_view w;
    MML_STRUCT_INIT(&w, mml_mesh_view);
    w.vertex_count = 4;
    w.triangle_count = 2;
    w.positions = wp;
    w.triangles = wf;

    mml_vec3 wn[4];
    mml_compute_vertex_normals(&w, NULL, wn, 4, NULL);
    mml_vec3 big = vec_cross(vec_sub(wp[1], wp[0]),
                             vec_sub(wp[3], wp[0]));
    mml_vec3 small = vec_cross(vec_sub(wp[1], wp[0]),
                               vec_sub(wp[2], wp[0]));
    check("bigger_dominates", vec_length(big) > vec_length(small));
    check("tilt_toward_big",
          fabsf(wn[0].y) > fabsf(wn[0].z) * 0.5f);

    mml_vec3 ep[1] = { { 0, 0, 0 } };
    mml_mesh_view empty;
    MML_STRUCT_INIT(&empty, mml_mesh_view);
    empty.vertex_count = 1;
    empty.positions = ep;
    mml_vec3 en[1];
    check("no_faces_rejected",
          mml_compute_vertex_normals(&empty, NULL, en, 1, NULL) < 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
