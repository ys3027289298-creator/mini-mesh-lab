/* tangent: reconstruct per-vertex tangent frames from UV gradients.
 *
 * Frame convention: T = dP/du direction, B = dP/dv direction, N = area-weighted
 * vertex normal.  A valid frame is orthonormal AND right-handed: cross(T,B)=+N.
 *
 * Raw UV gradients are only a direction source.  T is the tangential part of the
 * accumulated dP/du; B is then cross(N,T).  That makes both checks exact by
 * construction; the raw dP/dv is only used to detect a mirrored UV layout.
 *
 * When a vertex has no tangential dP/du (collapsed UV island, pole, all incident
 * UV triangles degenerate), the frame is approximated from a UV axis projected
 * into the tangent plane and marked TAN_FALLBACK.  A mirrored UV layout cannot
 * be right-handed while keeping T along dP/du: B is flipped to satisfy the
 * handedness requirement and the vertex is marked TAN_UV_MIRRORED.
 *
 * The build also runs a control mode (raw normalized gradients, no
 * orthonormalization) so the test can show what the fix buys and whether the
 * orthogonality residual tracks UV anisotropy.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 256
#define MAXF 512
#define PI_F 3.14159265f
#define RAD_TO_DEG 57.29578f

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

static V3 v3(float x, float y, float z) { V3 r = { x, y, z }; return r; }
static V3 vadd(V3 p, V3 q) { return v3(p.x+q.x, p.y+q.y, p.z+q.z); }
static V3 vsub(V3 p, V3 q) { return v3(p.x-q.x, p.y-q.y, p.z-q.z); }
static V3 vmul(V3 p, float s) { return v3(p.x*s, p.y*s, p.z*s); }
static V3 vcross(V3 p, V3 q) { return v3(p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x); }
static float vdot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float vlen(V3 p) { return sqrtf(vdot(p,p)); }
static V3 vnorm(V3 p) { float l = vlen(p); return l > 0.0f ? vmul(p, 1.0f/l) : v3(0,0,0); }

typedef struct {
    V3 v[MAXV]; float uv[MAXV][2];
    Tri f[MAXF]; int nv, nf;
} Mesh;

typedef struct {
    V3 t[MAXV], b[MAXV], n[MAXV];
    V3 uv_axis[MAXV];          /* area-weighted normalized dP/du, per vertex */
    unsigned flags[MAXV];
    int fallback;              /* vertices that used the UV-axis approximation */
    int mirrored;              /* vertices whose raw dP/dv opposes cross(N,T) */
    int degenerate_faces;      /* faces with zero UV area, not accumulated */
    float aniso_mean, aniso_max;   /* UV metric singular-value ratio per face */
    float shear_mean;              /* |b|/sqrt(ac): normalized UV cross-term */
} Frames;

enum { TAN_FALLBACK = 1u << 0, TAN_DEGENERATE_UV = 1u << 1, TAN_UV_MIRRORED = 1u << 2 };

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

static int mesh_add_vertex(Mesh* m, float x, float y, float z, float u, float v) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv] = v3(x, y, z);
    m->uv[m->nv][0] = u; m->uv[m->nv][1] = v;
    return m->nv++;
}

static int mesh_add_tri(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a < 0 || b < 0 || c < 0 || a >= m->nv || b >= m->nv || c >= m->nv) return -1;
    if (a == b || b == c || a == c) return -1;
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* Build frames.  raw=0: orthonormalized (production path);
 *             raw=1: raw normalized gradients (control, no orthonormalization).
 * Returns the number of fallback vertices, or -1 when the mesh is empty. */
static int tangent_build(const Mesh* m, int raw, Frames* fr) {
    if (m->nv <= 0 || m->nf <= 0) return -1;
    memset(fr, 0, sizeof(*fr));

    V3 g_u[MAXV], g_b[MAXV], f_n[MAXV];
    memset(g_u, 0, sizeof(g_u)); memset(g_b, 0, sizeof(g_b)); memset(f_n, 0, sizeof(f_n));
    float aniso_sum = 0.0f, shear_sum = 0.0f; int aniso_n = 0, shear_n = 0;

    for (int i = 0; i < m->nf; i++) {
        const Tri* f = &m->f[i];
        V3 p0 = m->v[f->a], p1 = m->v[f->b], p2 = m->v[f->c];
        V3 e1 = vsub(p1, p0), e2 = vsub(p2, p0);
        V3 nf = vcross(e1, e2);              /* length = 2*area */
        float area = 0.5f * vlen(nf);
        f_n[f->a] = vadd(f_n[f->a], nf);
        f_n[f->b] = vadd(f_n[f->b], nf);
        f_n[f->c] = vadd(f_n[f->c], nf);

        float d1u = m->uv[f->b][0]-m->uv[f->a][0], d1v = m->uv[f->b][1]-m->uv[f->a][1];
        float d2u = m->uv[f->c][0]-m->uv[f->a][0], d2v = m->uv[f->c][1]-m->uv[f->a][1];
        float det = d1u*d2v - d2u*d1v;
        if (fabsf(det) < 1e-12f) {
            fr->degenerate_faces++;
            fr->flags[f->a] |= TAN_DEGENERATE_UV;
            fr->flags[f->b] |= TAN_DEGENERATE_UV;
            fr->flags[f->c] |= TAN_DEGENERATE_UV;
            continue;
        }
        float inv = 1.0f / det;
        /* dP/du = (e1*d2v - e2*d1v)/det , dP/dv = (e2*d1u - e1*d2u)/det */
        V3 tu = vmul(vsub(vmul(e1, d2v), vmul(e2, d1v)), inv);
        V3 tv = vmul(vsub(vmul(e2, d1u), vmul(e1, d2u)), inv);
        V3 tuw = vmul(tu, area), tvw = vmul(tv, area);
        g_u[f->a] = vadd(g_u[f->a], tuw); g_u[f->b] = vadd(g_u[f->b], tuw); g_u[f->c] = vadd(g_u[f->c], tuw);
        g_b[f->a] = vadd(g_b[f->a], tvw); g_b[f->b] = vadd(g_b[f->b], tvw); g_b[f->c] = vadd(g_b[f->c], tvw);

        float tl = vlen(tu);
        if (tl > 1e-12f) {
            V3 au = vmul(tu, area / tl);
            fr->uv_axis[f->a] = vadd(fr->uv_axis[f->a], au);
            fr->uv_axis[f->b] = vadd(fr->uv_axis[f->b], au);
            fr->uv_axis[f->c] = vadd(fr->uv_axis[f->c], au);
        }
        /* UV anisotropy: sqrt of the eigenvalue ratio of the first fundamental
         * form induced by UV, i.e. how unevenly the mapping stretches the face. */
        float a = vdot(tu, tu), b2 = vdot(tu, tv), c = vdot(tv, tv);
        float disc = sqrtf((a-c)*(a-c) + 4.0f*b2*b2);
        float lmax = 0.5f*(a+c+disc), lmin = 0.5f*(a+c-disc);
        if (lmin > 1e-12f) {
            float ratio = sqrtf(lmax / lmin);
            aniso_sum += ratio; aniso_n++;
            if (ratio > fr->aniso_max) fr->aniso_max = ratio;
        }
        /* cos of the angle between the two raw UV gradients.  This, not the
         * anisotropy ratio, is what a non-orthonormalized reconstruction sees:
         * independently scaling dP/du or dP/dv (pure axis stretch) rescales a
         * and c but leaves b/sqrt(ac) unchanged. */
        float ac = sqrtf(a * c);
        if (ac > 1e-12f) { shear_sum += fabsf(b2) / ac; shear_n++; }
    }
    fr->aniso_mean = aniso_n > 0 ? aniso_sum / (float)aniso_n : 0.0f;
    fr->shear_mean = shear_n > 0 ? shear_sum / (float)shear_n : 0.0f;

    for (int i = 0; i < m->nv; i++) {
        V3 n = vnorm(f_n[i]);
        if (vlen(f_n[i]) <= 0.0f) n = v3(0, 0, 1);   /* isolated / zero-area fan */
        V3 gu = g_u[i], gb = g_b[i];
        V3 t = vsub(gu, vmul(n, vdot(n, gu)));
        V3 s = vsub(gb, vmul(n, vdot(n, gb)));
        float lt = vlen(t), ls = vlen(s);
        V3 T, B;

        if (lt < 1e-9f) {
            /* no tangential UV gradient: approximate from a UV axis, mark it */
            V3 axis = fabsf(n.x) < 0.9f ? v3(1,0,0) : v3(0,1,0);
            T = vnorm(vsub(axis, vmul(n, vdot(n, axis))));
            B = vcross(n, T);
            fr->fallback++;
            fr->flags[i] |= TAN_FALLBACK;
        } else {
            T = vmul(t, 1.0f / lt);
            if (raw) {
                B = ls > 1e-9f ? vmul(s, 1.0f / ls) : vcross(n, T);
            } else {
                B = vcross(n, T);                 /* exact T.B = 0 and cross(T,B)=+N */
                if (ls > 1e-9f && vdot(B, s) < 0.0f) {
                    fr->mirrored++;
                    fr->flags[i] |= TAN_UV_MIRRORED;   /* UV is mirrored; RHS wins */
                }
            }
        }
        fr->t[i] = T; fr->b[i] = B; fr->n[i] = n;
    }
    return fr->fallback;
}

/* residual: max over vertices of |T.B|, |T.N|, |B.N|, ||T|-1|, ||B|-1|.
 * rhs: fraction of vertices with cross(T,B).N > 0.
 * angle: mean angle between T and the vertex UV main axis (dP/du), degrees. */
static void frame_stats(const Frames* fr, int nv, float* out_res, float* out_rhs, float* out_ang) {
    float res = 0.0f, ang_sum = 0.0f;
    int rhs = 0, nang = 0;
    for (int i = 0; i < nv; i++) {
        V3 t = fr->t[i], b = fr->b[i], n = fr->n[i];
        float r = fabsf(vdot(t, b));
        if (fabsf(vdot(t, n)) > r) r = fabsf(vdot(t, n));
        if (fabsf(vdot(b, n)) > r) r = fabsf(vdot(b, n));
        if (fabsf(vlen(t) - 1.0f) > r) r = fabsf(vlen(t) - 1.0f);
        if (fabsf(vlen(b) - 1.0f) > r) r = fabsf(vlen(b) - 1.0f);
        if (r > res) res = r;
        if (vdot(vcross(t, b), n) > 0.0f) rhs++;

        V3 a = vsub(fr->uv_axis[i], vmul(n, vdot(n, fr->uv_axis[i])));
        float la = vlen(a), lt = vlen(t);
        if (la > 1e-9f && lt > 1e-9f) {
            float c = vdot(t, a) / (lt * la);
            if (c > 1.0f) c = 1.0f;
            if (c < -1.0f) c = -1.0f;
            ang_sum += acosf(c) * RAD_TO_DEG;
            nang++;
        }
    }
    *out_res = res;
    *out_rhs = nv > 0 ? (float)rhs / (float)nv : 0.0f;
    *out_ang = nang > 0 ? ang_sum / (float)nang : -1.0f;
}

/* ---- models: a bumpy grid with three different UV layouts ---- */
enum { UV_ISO, UV_STRETCH, UV_SHEAR_MIRROR };

static void build_grid(Mesh* m, int layout, int n, float amp) {
    mesh_init(m);
    for (int j = 0; j <= n; j++) {
        for (int i = 0; i <= n; i++) {
            float x = (float)i / (float)n, y = (float)j / (float)n;
            float u, v;
            switch (layout) {
            case UV_STRETCH:      u = 8.0f * x;            v = y; break;
            case UV_SHEAR_MIRROR: u = 1.0f - (x + 0.8f*y); v = y; break;
            default:              u = x;                   v = y; break;
            }
            float z = amp * sinf(PI_F * x) * sinf(PI_F * y);
            mesh_add_vertex(m, x, y, z, u, v);
        }
    }
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            int a = j*(n+1) + i, b = a + 1, c = a + n + 2, d = a + n + 1;
            mesh_add_tri(m, a, b, c);
            mesh_add_tri(m, a, c, d);
        }
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    enum { NGRID = 12, NMODEL = 3 };
    const char* names[NMODEL] = { "iso", "stretch8", "shear_mirror" };
    const int layouts[NMODEL] = { UV_ISO, UV_STRETCH, UV_SHEAR_MIRROR };

    Mesh mesh[NMODEL]; Frames fx[NMODEL], ctl[NMODEL];
    float res[NMODEL], rhs[NMODEL], ang[NMODEL];
    float c_res[NMODEL], c_rhs[NMODEL], c_ang[NMODEL];

    int built = 1;
    for (int k = 0; k < NMODEL; k++) {
        build_grid(&mesh[k], layouts[k], NGRID, 0.25f);
        if (tangent_build(&mesh[k], 0, &fx[k]) < 0) built = 0;
        if (tangent_build(&mesh[k], 1, &ctl[k]) < 0) built = 0;
        frame_stats(&fx[k], mesh[k].nv, &res[k], &rhs[k], &ang[k]);
        frame_stats(&ctl[k], mesh[k].nv, &c_res[k], &c_rhs[k], &c_ang[k]);
    }
    check("built_three_models", built);

    /* --- the nine numbers: orthogonality residual, right-handed fraction,
     *     angle to the reconstructed tangent vs the UV main axis --- */
    printf("\nreconstructed model         orth_res  rhs_frac  angle_deg\n");
    for (int k = 0; k < NMODEL; k++)
        printf("%-19s %.3e  %.4f  %.4f\n", names[k], res[k], rhs[k], ang[k]);

    check("orth_res_iso", res[0] < 1e-5f);
    check("orth_res_stretch", res[1] < 1e-5f);
    check("orth_res_shear_mirror", res[2] < 1e-5f);
    check("rhs_all", rhs[0] == 1.0f && rhs[1] == 1.0f && rhs[2] == 1.0f);
    check("angle_all", ang[0] < 0.5f && ang[1] < 0.5f && ang[2] < 0.5f);
    check("no_fallback_on_clean_meshes",
          fx[0].fallback == 0 && fx[1].fallback == 0 && fx[2].fallback == 0);
    check("mirror_flagged", fx[2].mirrored == mesh[2].nv &&
                            (fx[2].flags[0] & TAN_UV_MIRRORED) != 0);

    /* --- control: raw normalized UV gradients, no orthonormalization.
     *     Residual and handedness are what a naive reconstruction would give. --- */
    printf("\ncontrol: raw UV gradients, no orthonormalization\n");
    printf("model         aniso_mean aniso_max  shear  orth_res  rhs_frac  angle_deg\n");
    for (int k = 0; k < NMODEL; k++)
        printf("%-13s %8.3f  %8.3f  %6.3f  %.3e  %.4f  %.4f\n",
               names[k], ctl[k].aniso_mean, ctl[k].aniso_max, ctl[k].shear_mean,
               c_res[k], c_rhs[k], c_ang[k]);

    int max_aniso = 0, max_raw_res = 0;
    for (int k = 1; k < NMODEL; k++) {
        if (ctl[k].aniso_mean > ctl[max_aniso].aniso_mean) max_aniso = k;
        if (c_res[k] > c_res[max_raw_res]) max_raw_res = k;
    }
    printf("\nmax_aniso_model=%s max_raw_res_model=%s stretch_has_max_raw_res=%d\n",
           names[max_aniso], names[max_raw_res],
           max_raw_res == 1 ? 1 : 0);
    printf("aniso_vs_residual: %s aniso=%.3f raw_res=%.3e | %s aniso=%.3f raw_res=%.3e | %s aniso=%.3f raw_res=%.3e\n",
           names[0], ctl[0].aniso_mean, c_res[0],
           names[1], ctl[1].aniso_mean, c_res[1],
           names[2], ctl[2].aniso_mean, c_res[2]);

    /* Same texture across models: where the UV main axis is the same (iso and
     * axis-stretch both map u along +x), the reconstructed T must be the same
     * direction despite the 8x UV scale difference. */
    float cross_max = 0.0f;
    for (int i = 0; i < mesh[0].nv; i++) {
        float c = vdot(fx[0].t[i], fx[1].t[i]);
        if (c > 1.0f) c = 1.0f;
        if (c < -1.0f) c = -1.0f;
        float d = acosf(c) * RAD_TO_DEG;
        if (d > cross_max) cross_max = d;
    }
    printf("cross_model_tangent_angle_iso_vs_stretch=%.4f deg\n", cross_max);
    check("same_tangent_across_uv_scale", cross_max < 0.5f);

    /* Axis-aligned stretch only rescales the two gradients; normalization erases
     * the scale, so the orthogonality residual cannot grow with anisotropy.
     * Residual is driven by the shear term, not by the anisotropy ratio. */
    check("raw_aniso_max_is_stretch", max_aniso == 1);
    check("raw_stretch_not_max_residual", max_raw_res == 2);
    check("raw_axis_stretch_leaves_residual", fabsf(c_res[0] - c_res[1]) < 1e-3f);
    check("raw_shear_dominates_residual", c_res[2] > 5.0f * c_res[0]);
    check("raw_mirror_breaks_handedness", c_rhs[2] < 0.5f);
    /* the residual tracks the shear term, not the anisotropy ratio */
    check("residual_tracks_shear_not_aniso",
          fabsf(c_res[2] - ctl[2].shear_mean) < 0.1f &&
          ctl[0].aniso_mean < 2.0f && ctl[1].aniso_mean > 7.0f &&
          fabsf(c_res[0] - c_res[1]) < 1e-3f);

    /* --- fallback: collapsed UV at one vertex must degrade and be marked --- */
    Mesh d; mesh_init(&d);
    mesh_add_vertex(&d, 0,0,0, 0,0);
    mesh_add_vertex(&d, 1,0,0, 1,0);
    mesh_add_vertex(&d, 1,1,0, 1,1);
    mesh_add_vertex(&d, 0,1,0, 0,1);
    mesh_add_vertex(&d, 0.5f,0.5f,0.5f, 0,0);   /* UV collapsed onto vertex 0 */
    mesh_add_tri(&d, 0,1,2); mesh_add_tri(&d, 0,2,3);
    mesh_add_tri(&d, 0,1,4);                     /* zero UV area */
    Frames df;
    int fb = tangent_build(&d, 0, &df);
    check("fallback_count", fb == 1);
    check("fallback_flagged", (df.flags[4] & TAN_FALLBACK) != 0);
    check("degenerate_uv_flagged", (df.flags[4] & TAN_DEGENERATE_UV) != 0);
    float dr, dh, da;
    frame_stats(&df, d.nv, &dr, &dh, &da);
    check("fallback_residual", dr < 1e-5f);
    check("fallback_rhs", dh == 1.0f);
    check("fallback_t_in_plane", fabsf(vdot(df.t[4], df.n[4])) < 1e-6f);

    Mesh e; mesh_init(&e);
    mesh_add_vertex(&e, 0,0,0, 0,0);
    Frames ef;
    check("empty_rejected", tangent_build(&e, 0, &ef) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
