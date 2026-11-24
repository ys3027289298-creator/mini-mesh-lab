
/* skin: linear-blend skinning, weight normalization, and elbow shrink probe.
 *
 * Lab fixture: a two-segment arm. Bone 0 is the upper arm (identity in every
 * pose); bones 1..7 are the forearm family and share one rotation Rz(theta).
 * The elbow ring is weighted on all 7 forearm bones with a raw sum s < 1, so
 * the classic real-time LBS formula q = sum(w_i * T_i * p) (no homogeneous
 * divide) pulls the ring toward the origin. See README for the derivation.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846f
#define NB     8    /* bone 0 upper arm; bones 1..7 forearm family */
#define NR     12   /* vertices on the elbow cross-section ring */
#define RADIUS 0.5f
#define LEN    1.0f /* upper arm / forearm length, joint at x = LEN */
#define CAP    4    /* hardware influence limit under test */
#define W_MAIN 0.80f
#define W_AUX  0.02f /* six auxiliary bones: 6 * 0.02 = 0.12 */
#define S_RAW  (W_MAIN + 6 * W_AUX)            /* 0.92 */
#define S_KEPT (W_MAIN + (CAP - 1) * W_AUX)    /* 0.86 after top-4 cap */

typedef struct { float x, y, z; } V3;

enum { MODE_RAW = 0, MODE_NORM, MODE_TRUNC, MODE_TRUNC_RN, NMODE };
static const char* mode_name[NMODE] = { "raw", "normalized", "trunc4", "trunc4_renorm" };

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* bone transforms: upper arm identity, forearm family Rz(theta) about origin */
static void bone_pose(int bone, float th, const V3* p, V3* q) {
    if (bone == 0) { *q = *p; return; }
    float c = cosf(th), s = sinf(th);
    q->x = c * p->x - s * p->y;
    q->y = s * p->x + c * p->y;
    q->z = p->z;
}

static int pick_top(const float* w, int cap, int* idx) {
    float tmp[NB];
    memcpy(tmp, w, sizeof(tmp));
    int k = 0;
    for (int c = 0; c < cap; c++) {
        int best = -1;
        for (int i = 0; i < NB; i++)
            if (tmp[i] > 0.0f && (best < 0 || tmp[i] > tmp[best])) best = i;
        if (best < 0) break;
        idx[k++] = best;
        tmp[best] = -1.0f;
    }
    return k;
}

/* LBS one vertex. Reports the weight sum actually used and influence count. */
static void skin_vertex(const float* w0, const V3* p, float th, int mode,
                        V3* q, float* sum_out, int* inf_out) {
    float w[NB];
    memcpy(w, w0, sizeof(w));
    int inf = 0;
    float s = 0.0f;
    for (int i = 0; i < NB; i++) if (w[i] > 0.0f) { inf++; s += w[i]; }

    if (mode == MODE_NORM) {
        for (int i = 0; i < NB; i++) w[i] /= s;
        s = 1.0f;
    } else if (mode == MODE_TRUNC || mode == MODE_TRUNC_RN) {
        int idx[NB];
        int k = pick_top(w0, CAP, idx);
        float kept = 0.0f;
        for (int c = 0; c < k; c++) kept += w0[idx[c]];
        memset(w, 0, sizeof(w));
        for (int c = 0; c < k; c++)
            w[idx[c]] = (mode == MODE_TRUNC_RN) ? w0[idx[c]] / kept : w0[idx[c]];
        inf = k;
        s = (mode == MODE_TRUNC_RN) ? 1.0f : kept;
    }

    float ax = 0.0f, ay = 0.0f, az = 0.0f;
    for (int i = 0; i < NB; i++) {
        if (w[i] <= 0.0f) continue;
        V3 t;
        bone_pose(i, th, p, &t);
        ax += w[i] * t.x; ay += w[i] * t.y; az += w[i] * t.z;
    }
    q->x = ax; q->y = ay; q->z = az;
    *sum_out = s;
    *inf_out = inf;
}

typedef struct { float perim, avg_d; } Meas;

/* elbow ring: axis is the posed forearm line through the elbow joint e */
static Meas measure(const V3 ring[NR], float th) {
    float ax = cosf(th), ay = sinf(th);
    float ex = LEN * ax, ey = LEN * ay;
    float perim = 0.0f, dsum = 0.0f;
    for (int j = 0; j < NR; j++) {
        const V3* a = &ring[j];
        const V3* b = &ring[(j + 1) % NR];
        float dx = a->x - b->x, dy = a->y - b->y, dz = a->z - b->z;
        perim += sqrtf(dx * dx + dy * dy + dz * dz);

        float vx = a->x - ex, vy = a->y - ey, vz = a->z;
        float dot = vx * ax + vy * ay;
        float px = vx - dot * ax, py = vy - dot * ay, pz = vz;
        dsum += sqrtf(px * px + py * py + pz * pz);
    }
    Meas m;
    m.perim = perim;
    m.avg_d = dsum / NR;
    return m;
}

int main(void) {
    V3 elbow_rest[NR], elbow_sk[NMODE][NR];
    float w_elbow[NR][NB];

    memset(w_elbow, 0, sizeof(w_elbow));
    for (int j = 0; j < NR; j++) {
        float phi = 2.0f * PI * j / NR;
        elbow_rest[j] = (V3){ LEN, RADIUS * cosf(phi), RADIUS * sinf(phi) };
        w_elbow[j][1] = W_MAIN;
        for (int b = 2; b < NB; b++) w_elbow[j][b] = W_AUX;
    }

    /* rigid reference: polygon perimeter and radius are pose-independent */
    float p0 = 2.0f * NR * RADIUS * sinf(PI / NR);
    const float deg[3] = { 0.0f, 45.0f, 90.0f };
    const float th[3]  = { 0.0f, PI / 4.0f, PI / 2.0f };

    check("raw_sum_below_one", fabsf(S_RAW - 0.92f) < 1e-6f);
    check("raw_influences_over_cap", 7 > CAP);
    check("kept_sum_after_cap", fabsf(S_KEPT - 0.86f) < 1e-6f);

    float ratio_p[NMODE][3], ratio_d[NMODE][3];
    float used_sum[NMODE]; int used_inf[NMODE];

    printf("\n%-14s %5s %10s %10s %9s %9s %7s\n",
           "mode", "angle", "perim", "avg_d", "P/P0", "d/r", "wsum");
    for (int md = 0; md < NMODE; md++) {
        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < NR; j++)
                skin_vertex(w_elbow[j], &elbow_rest[j], th[k], md,
                            &elbow_sk[md][j], &used_sum[md], &used_inf[md]);
            Meas m = measure(elbow_sk[md], th[k]);
            ratio_p[md][k] = m.perim / p0;
            ratio_d[md][k] = m.avg_d / RADIUS;
            printf("%-14s %5.0f %10.6f %10.6f %9.5f %9.5f %7.4f\n",
                   mode_name[md], deg[k], m.perim, m.avg_d,
                   ratio_p[md][k], ratio_d[md][k], used_sum[md]);
        }
    }

    /* relation q_raw = s * q_norm: every skinned length scales by the sum */
    float max_err = 0.0f;
    for (int k = 0; k < 3; k++)
        for (int j = 0; j < NR; j++) {
            V3* a = &elbow_sk[MODE_RAW][j];
            V3* b = &elbow_sk[MODE_NORM][j];
            float ex = a->x - S_RAW * b->x, ey = a->y - S_RAW * b->y,
                  ez = a->z - S_RAW * b->z;
            float e = sqrtf(ex * ex + ey * ey + ez * ez);
            if (e > max_err) max_err = e;
        }
    printf("\nrelation_q_eq_s_qnorm_maxerr=%.2e\n", max_err);
    check("relation_q_eq_s_qnorm", max_err < 1e-5f);

    /* raw: shrink factor equals the weight sum at every pose */
    for (int k = 0; k < 3; k++) {
        check("raw_scales_by_sum",
              fabsf(ratio_p[MODE_RAW][k] - S_RAW) < 1e-4f &&
              fabsf(ratio_d[MODE_RAW][k] - S_RAW) < 1e-4f);
    }
    /* normalized: 90-degree perimeter must match the rest pose */
    check("norm_no_shrink_at_90",
          fabsf(ratio_p[MODE_NORM][2] - 1.0f) < 1e-5f &&
          fabsf(ratio_d[MODE_NORM][2] - 1.0f) < 1e-5f);
    check("norm_sums_to_one", fabsf(used_sum[MODE_NORM] - 1.0f) < 1e-6f);

    /* truncation without renormalization: the failure path comes back */
    check("trunc_count_is_cap", used_inf[MODE_TRUNC] == CAP);
    check("trunc_shrinks_at_90",
          ratio_p[MODE_TRUNC][2] < 0.95f &&
          fabsf(ratio_p[MODE_TRUNC][2] - S_KEPT) < 1e-4f);
    check("trunc_renorm_repairs",
          fabsf(ratio_p[MODE_TRUNC_RN][2] - 1.0f) < 1e-5f &&
          fabsf(used_sum[MODE_TRUNC_RN] - 1.0f) < 1e-6f);

    printf("\n--- truncation diagnostic at the elbow ---\n");
    printf("actual_influences=%d cap=%d kept_influences=%d\n", 7, CAP, CAP);
    printf("raw_sum=%.4f normalized_sum=%.4f\n", S_RAW, 1.0f);
    printf("trunc_kept_sum_no_renorm=%.4f trunc_kept_sum_renorm=%.4f\n",
           S_KEPT, 1.0f);
    printf("trunc90_perim_vs_rest=%.4f (明显小于1: %s)\n",
           ratio_p[MODE_TRUNC][2], ratio_p[MODE_TRUNC][2] < 0.95f ? "yes" : "no");

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
