
/* seam: proxy metrics for "no visible seam at the bevel".
 * Textually includes the production normal code so the area-weighted
 * vertex_normals() under test is exactly normal/src/normal.c.
 * Compile with -Dmain=normal_demo_main to drop that file's main(). */
#include "../../normal/src/normal.c"
/* keep this TU's entry point named main despite -Dmain=normal_demo_main */
#undef main

#include <stdlib.h>

#define RADIUS 1.0f
#define FLAT_LEN 2.0f
#define EXTRUDE 1.0f
#define M_SAMPLES 400          /* fine samples along the profile line */
#define JND 0.01f              /* ~1 Weber fraction, full-scale luminance */

/* ---- profile: flat1 (-x) -> quarter-circle bevel -> flat2 (+y) -------- */
typedef struct { V3 p; float s; int is_arc; } PNode;

static V3 v2(float x, float y) { V3 r = { x, y, 0.0f }; return r; }

/* N bevel segments, uniform angle if crowd==0, quadratic crowding if crowd==1.
 * Fillet: circle centre (0,1), r=1, tangent to flat1 y=0 at (0,0) and
 * flat2 x=1 at (1,1). Nodes: [(-2,0)] [(0,0) seam] arc j=1..N [(1,3)] */
static int build_profile(PNode* nodes, int N, int crowd) {
    int k = 0;
    nodes[k].p = v2(-FLAT_LEN, 0); nodes[k].s = 0; nodes[k].is_arc = 0; k++;
    nodes[k].p = v2(0.0f, 0.0f);   nodes[k].s = FLAT_LEN; nodes[k].is_arc = 0; k++;
    for (int j = 1; j <= N; j++) {
        float t = (float)j / (float)N;
        float u = crowd ? t * t : t;
        float phi = (-0.5f + 0.5f * u) * 3.14159265f;   /* -pi/2 .. 0 */
        PNode nd;
        nd.p = v2(RADIUS * cosf(phi), 1.0f + RADIUS * sinf(phi));
        nd.s = nodes[k - 1].s + len(sub(nd.p, nodes[k - 1].p));
        nd.is_arc = (j < N);      /* j==N is the flat2 boundary node */
        nodes[k++] = nd;
    }
    nodes[k].p = v2(1.0f, 1.0f + FLAT_LEN);
    nodes[k].s = nodes[k - 1].s + FLAT_LEN;
    nodes[k].is_arc = 0; k++;
    return k;
}

/* extrude profile along z, two triangles per profile segment */
static int build_mesh(const PNode* nodes, int np, Mesh* m) {
    memset(m, 0, sizeof(*m));
    int base[MAXV];
    for (int i = 0; i < np; i++) {
        int a = add_v(m, nodes[i].p.x, nodes[i].p.y, 0.0f);
        int b = add_v(m, nodes[i].p.x, nodes[i].p.y, EXTRUDE);
        if (a < 0 || b < 0) return -1;
        base[i] = a;
    }
    for (int i = 0; i + 1 < np; i++) {
        int a = base[i], b = base[i + 1], c = base[i + 1] + 1, d = base[i] + 1;
        if (add_f(m, a, b, c) < 0 || add_f(m, a, c, d) < 0) return -1;
    }
    return 0;
}

/* ---- per-node normals for the three schemes --------------------------- */
static void node_normals_flat(const PNode* nodes, int np, V3* nn) {
    for (int i = 0; i + 1 < np; i++) {
        V3 t = sub(nodes[i + 1].p, nodes[i].p);
        float l = len(t);
        /* outward 2D normal of tangent (t_x,t_y): (t_y,-t_x)/l */
        nn[i].x = t.y / l; nn[i].y = -t.x / l; nn[i].z = 0;
    }
    nn[np - 1] = nn[np - 2];
}

static void node_normals_area(const Mesh* m, int np, V3* nn) {
    V3 vn[MAXV];
    vertex_normals(m, vn);                 /* production call */
    for (int i = 0; i < np; i++) {         /* average the two extrusion verts */
        V3 a = vn[2 * i], b = vn[2 * i + 1];
        nn[i] = norm(v2(a.x + b.x, a.y + b.y));
        nn[i].z = 0;
        nn[i] = norm(nn[i]);
    }
}

/* reference field and fix coincide: analytic radial normals on arc nodes,
 * face normals on flats. Boundary nodes satisfy both by construction. */
static void node_normals_radial(const PNode* nodes, int np, V3* nn) {
    node_normals_flat(nodes, np, nn);
    for (int i = 0; i < np; i++) {
        if (nodes[i].is_arc) {
            V3 rad = sub(nodes[i].p, v2(0, 1.0f));
            nn[i] = norm(rad);
        }
    }
}

/* ---- metric A: normal-gradient jump vs designed field ----------------- */
/* g_i = |(N_{i+1}-N_i)/h_i - (N_i-N_{i-1})/h_{i-1}|, h in units of r.
 * R_i = max(0, g_i - g_ref_i). Also report raw g to show why the naive
 * second difference alone cannot gate a real fillet. */
static float grad_jump(const PNode* nd, const V3* nn, int i) {
    float h1 = (nd[i].s - nd[i - 1].s) / RADIUS;
    float h2 = (nd[i + 1].s - nd[i].s) / RADIUS;
    V3 gl = sub(nn[i], nn[i - 1]); gl.x /= h1; gl.y /= h1; gl.z /= h1;
    V3 gr = sub(nn[i + 1], nn[i]); gr.x /= h2; gr.y /= h2; gr.z /= h2;
    return len(sub(gr, gl));
}

static void metric_a(const PNode* nd, int np, const V3* nn, const V3* ref,
                     int seam_lo, int seam_hi, int flat_mode,
                     float* ngj_seam, float* ngj_max, float* raw_seam,
                     float* hard_gap, float* dev_rms) {
    float gj_s = 0, gj_m = 0, raw_s = 0, gap = 0, acc = 0;
    for (int i = 1; i + 1 < np; i++) {
        float g = grad_jump(nd, nn, i);
        float gr = grad_jump(nd, ref, i);
        float res = g - gr; if (res < 0) res = 0;
        if (res > gj_m) gj_m = res;
        if (i == seam_lo || i == seam_hi) {
            if (res > gj_s) gj_s = res;
            if (g > raw_s) raw_s = g;
        }
        V3 d = sub(nn[i], ref[i]);
        acc += dot(d, d);
    }
    if (flat_mode) {
        gap = len(sub(nn[seam_lo], nn[seam_lo - 1]));
        float g2 = len(sub(nn[seam_hi], nn[seam_hi - 1]));
        if (g2 > gap) gap = g2;
    }
    *ngj_seam = gj_s; *ngj_max = gj_m; *raw_seam = raw_s;
    *hard_gap = gap; *dev_rms = sqrtf(acc / (np - 2));
}

/* interpolated unit normal at fine arc-length position s */
static V3 normal_at(const PNode* nd, const V3* nn, int np, float s, int flat_mode) {
    int i = 0;
    while (i + 2 < np && nd[i + 1].s < s) i++;
    float h = nd[i + 1].s - nd[i].s;
    float t = (s - nd[i].s) / h;
    if (flat_mode) return nn[i];
    V3 a = nn[i], b = nn[i + 1];
    V3 r = { a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, 0 };
    return norm(r);
}

/* ---- metric B: luminance corner energy along the sampling line --------
 * L_i = max(0,N_i.l); segment slope; J_i = |slope_i - slope_{i-1}|;
 * E_i = J_i * min(h): full-scale luminance fraction over one segment.
 * Raw E fires on every facet of a tessellated smooth fillet, so the gated
 * quantity is the RESIDUAL E_i - E_ref_i against the designed field on the
 * same tessellation. B0 = C0 luminance jump (flat shading only). */
static void corner_energy(const PNode* nd, int np, const V3* nn, V3 l,
                          int flat_mode, float* E, float* b0) {
    float L[MAXV], slope[MAXV];
    for (int i = 0; i < np; i++) {
        float d = dot(nn[i], l);
        L[i] = d > 0 ? d : 0;
    }
    for (int i = 0; i + 1 < np; i++)
        slope[i] = (L[i + 1] - L[i]) / ((nd[i + 1].s - nd[i].s) / RADIUS);
    float jump_max = 0;
    E[0] = 0;
    for (int i = 1; i + 1 < np; i++) {
        float w = fminf(nd[i].s - nd[i - 1].s, nd[i + 1].s - nd[i].s) / RADIUS;
        E[i] = fabsf(slope[i] - slope[i - 1]) * w;
        if (flat_mode) {
            float jmp = fabsf(L[i] - L[i - 1]);
            if (jmp > jump_max) jump_max = jmp;
        }
    }
    E[np - 1] = 0;
    *b0 = jump_max;
}

static void metric_b(const PNode* nd, int np, const V3* nn, const V3* ref,
                     V3 l, int flat_mode,
                     float* lce_res, int* corners, float* b0) {
    float E[MAXV], Er[MAXV], dummy;
    corner_energy(nd, np, nn, l, flat_mode, E, b0);
    corner_energy(nd, np, ref, l, 0, Er, &dummy);
    float e_max = 0; int c = 0;
    for (int i = 1; i + 1 < np; i++) {
        float r = E[i] - Er[i];
        if (r > e_max) e_max = r;
        if (r > JND) c++;
    }
    *lce_res = e_max; *corners = c;
}

static void run_case(const char* tag, int N, int crowd, FILE* csv) {
    PNode nd[MAXV];
    int np = build_profile(nd, N, crowd);
    Mesh m;
    build_mesh(nd, np, &m);

    V3 nf[MAXV], na[MAXV], nx[MAXV], ref[MAXV];
    node_normals_flat(nd, np, nf);
    node_normals_area(&m, np, na);
    node_normals_radial(nd, np, nx);
    node_normals_radial(nd, np, ref);

    int seam_lo = 1, seam_hi = np - 2;   /* bevel boundary profile nodes */

    V3 la = norm(v2(0.6f, -0.8f));       /* lights the whole 90deg span */
    const char* scheme[3] = { "flat", "area(current)", "fixed" };
    V3* sets[3] = { nf, na, nx };
    int isflat[3] = { 1, 0, 0 };

    printf("== case %s (N=%d%s) total_length=%.3f seam_nodes=%d,%d ==\n",
           tag, N, crowd ? ", non-uniform" : ", uniform", nd[np - 1].s, seam_lo, seam_hi);
    printf("%-15s %8s %8s %9s %9s %8s %7s %8s\n",
           "scheme", "NGJres", "NGJmax", "rawD2seam", "hardgap", "devRMS",
           "corners", "LCEres");
    for (int k = 0; k < 3; k++) {
        float gj_s, gj_m, raw, gap, rms;
        metric_a(nd, np, sets[k], ref, seam_lo, seam_hi, isflat[k],
                 &gj_s, &gj_m, &raw, &gap, &rms);
        float lce, b0; int corners;
        metric_b(nd, np, sets[k], ref, la, isflat[k], &lce, &corners, &b0);
        printf("%-15s %8.4f %8.4f %9.4f %9.4f %8.4f %7d %8.4f",
               scheme[k], gj_s, gj_m, raw, gap, rms, corners, lce);
        if (isflat[k]) printf("  B0_jump=%.4f", b0);
        printf("\n");
    }

    /* B blind-spot demo: rotate the light to grazing incidence on one seam.
     * The Lambert clamp zeroes that seam's luminance response, so LCEres
     * stops reporting the defect there while A (light-independent) stays
     * exactly as bad. This is why A, not B, is the gate. */
    V3 lg = norm(v2(0.0f, -1.0f));       /* parallel to flat1 normal (0,-1) */
    float ce_n[MAXV], ce_r[MAXV], dummy2;
    corner_energy(nd, np, na, lg, 0, ce_n, &dummy2);
    corner_energy(nd, np, ref, lg, 0, ce_r, &dummy2);
    float res_grazing = fmaxf(ce_n[seam_lo] - ce_r[seam_lo],
                              ce_n[seam_hi] - ce_r[seam_hi]);
    float ce_n2[MAXV], ce_r2[MAXV], dummy3;
    corner_energy(nd, np, na, la, 0, ce_n2, &dummy3);
    corner_energy(nd, np, ref, la, 0, ce_r2, &dummy3);
    float res_chosen = fmaxf(ce_n2[seam_lo] - ce_r2[seam_lo],
                             ce_n2[seam_hi] - ce_r2[seam_hi]);
    float gj_area, gjm, raw_a, gap_a, rms_a;
    metric_a(nd, np, na, ref, seam_lo, seam_hi, 0,
             &gj_area, &gjm, &raw_a, &gap_a, &rms_a);
    printf("blind_spot area(current): seam LCEres %.4f (normal light) -> %.4f "
           "(grazing light); NGJres stays %.4f either way\n\n",
           res_chosen, res_grazing, gj_area);

    /* sanity dump: normals actually present at the two seam nodes */
    printf("seam node normals (flat / area / reference):\n");
    for (int i = 0; i < np; i++) {
        if (i == seam_lo || i == seam_hi) {
            printf("  node %d: (%.4f,%.4f) (%.4f,%.4f) (%.4f,%.4f)\n", i,
                   nf[i].x, nf[i].y, na[i].x, na[i].y, ref[i].x, ref[i].y);
        }
    }
    printf("\n");

    float S = nd[np - 1].s;
    for (int q = 0; q <= M_SAMPLES; q++) {
        float s = S * q / M_SAMPLES;
        V3 vf = normal_at(nd, nf, np, s, 1);
        V3 va = normal_at(nd, na, np, s, 0);
        V3 vx = normal_at(nd, nx, np, s, 0);
        fprintf(csv, "%s,%.5f,%.5f,%.5f,%.5f\n", tag, s,
                fmaxf(0, dot(vf, la)), fmaxf(0, dot(va, la)), fmaxf(0, dot(vx, la)));
    }
}

int main(void) {
    FILE* csv = fopen("out/lum.csv", "w");
    if (!csv) { perror("out/lum.csv"); return 1; }
    fprintf(csv, "case,s,L_flat,L_area,L_fixed\n");
    run_case("uniform8", 8, 0, csv);
    run_case("nonuniform8", 8, 1, csv);
    fclose(csv);
    return 0;
}
