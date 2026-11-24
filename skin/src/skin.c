
/* skin: limit vertex influences to four bones, error-aware ordering.
 *
 * A vertex whose influence set is over the 4-bone budget is truncated:
 *   - candidates are sorted by  score = w_n * (r / R)^beta
 *       w_n : normalized bone weight  (weight factor)
 *       r   : distance from the vertex to that bone's joint (pivot)
 *       R   : model radius, so r/R is dimensionless
 *   - the top four survive, the rest are dropped
 *   - the kept weights are renormalized to sum to one
 *   - L = dropped raw weight / total raw weight is the per-vertex weight loss
 *
 * beta = 0 is plain weight-order pruning.  A weight-only key treats a bone
 * that barely touches the vertex but has a long lever as cheap to drop, even
 * though dropping it moves the vertex; beta > 0 prices that lever in.  beta
 * is not a constant of the algorithm: it is solved per model by minimizing
 * the measured truncation error (see solve_beta).
 *
 * Test rig: an elbow.  A cylinder of radius R is skinned to two axial bones
 * plus a ring of auxiliary "collar" bones around the joint, so vertices near
 * the elbow see well over four bones.  Bone density varies per model while
 * the surface tessellation is held fixed, so only the bone count changes.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 512
#define MAXB 16
#define MAXK 4

#define R_MESH   0.35f     /* cylinder radius */
#define RING_A   14        /* angular samples, fixed across models */
#define RING_X   20        /* axial samples, fixed across models */
#define X0       0.35f
#define X1       1.65f
#define JOINT_X  1.0f      /* elbow */
#define BEND     0.6f      /* radians the lower arm rotates */
#define WEIGHT_P 3.0f
#define WEIGHT_D 0.02f     /* softening so weight stays finite on the axis */
#define CUTOFF   0.50f     /* a bone influences a vertex only inside this */
#define COLLAR_R 0.15f     /* collar bones sit this far off the axis */
#define COLLAR_X0 0.55f
#define COLLAR_X1 1.45f

#define BETA_MAX 6.0f
#define BETA_STEP 0.01f

typedef struct { float x, y, z; } V3;
typedef struct { V3 p0, p1; float angle; } Bone;   /* rotate about p0, rz */

static V3 v3(float x, float y, float z) { V3 r = { x, y, z }; return r; }
static V3 vsub(V3 a, V3 b) { return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static V3 vadd(V3 a, V3 b) { return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static V3 vmul(V3 a, float s) { return v3(a.x*s, a.y*s, a.z*s); }
static float vdot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static float vlen(V3 a) { return sqrtf(vdot(a, a)); }
static V3 rotz(V3 p, V3 pivot, float a) {
    V3 d = vsub(p, pivot);
    float c = cosf(a), s = sinf(a);
    return v3(pivot.x + c*d.x - s*d.y, pivot.y + s*d.x + c*d.y, pivot.z + d.z);
}

/* distance from p to the bone segment */
static float seg_dist(V3 p, const Bone* b) {
    V3 d = vsub(b->p1, b->p0);
    float dd = vdot(d, d);
    float t = dd > 0.0f ? vdot(vsub(p, b->p0), d) / dd : 0.0f;
    if (t < 0.0f) t = 0.0f; else if (t > 1.0f) t = 1.0f;
    return vlen(vsub(p, v3(b->p0.x + d.x*t, b->p0.y + d.y*t, b->p0.z + d.z*t)));
}

typedef struct {
    const char* name;
    int density;
    int nv;
    V3 rest[MAXV];
    int nb;
    Bone bones[MAXB];
    float radius;
} Model;

static void build_model(Model* m, const char* name, int density) {
    memset(m, 0, sizeof(*m));
    m->name = name;
    m->density = density;

    /* surface: fixed tessellation regardless of bone density */
    for (int ix = 0; ix < RING_X; ix++) {
        float x = X0 + (X1 - X0) * (float)ix / (float)(RING_X - 1);
        for (int ia = 0; ia < RING_A; ia++) {
            float a = 6.2831853f * (float)ia / (float)RING_A;
            m->rest[m->nv++] = v3(x, R_MESH * cosf(a), R_MESH * sinf(a));
        }
    }

    /* two axial bones: upper arm static, lower arm bends at the elbow */
    m->bones[0] = (Bone){ v3(0,0,0), v3(JOINT_X,0,0), 0.0f };
    m->bones[1] = (Bone){ v3(JOINT_X,0,0), v3(2*JOINT_X,0,0), BEND };

    /* collar bones: parallel to the arm, offset, own pivots along the joint.
     * They swing hard against the elbow bend, so dropping one is a real
     * displacement and the ordering key actually matters. */
    int n_collar = density - 2;
    for (int i = 0; i < n_collar; i++) {
        float th = 6.2831853f * (float)i / (float)n_collar;
        float cy = COLLAR_R * cosf(th), cz = COLLAR_R * sinf(th);
        int seg = i & 1;    /* alternate which half of the arm the bone spans */
        float x0 = seg ? JOINT_X : COLLAR_X0;
        float x1 = seg ? COLLAR_X1 : JOINT_X;
        float ang = BEND * ((i & 1) ? -1.6f : 1.6f);
        m->bones[2 + i] = (Bone){ v3(x0, cy, cz), v3(x1, cy, cz), ang };
    }
    m->nb = 2 + n_collar;

    m->radius = 1e-6f;
    for (int i = 0; i < m->nv; i++) {
        float r = vlen(m->rest[i]);
        if (r > m->radius) m->radius = r;
    }
}

typedef struct {
    float w[MAXB];          /* final normalized weights, 0 where dropped */
    int   kept[MAXK];       /* kept bone indices */
    int   nkept;
    int   ninfl;            /* raw influence count */
    float loss;             /* dropped / total, in [0,1] */
    float score[MAXB];      /* ordering key per bone, 0 where unused */
} Skin;

/* raw inverse-distance weights; returns the raw influence count */
static int raw_weights(const Model* m, V3 p, float* w) {
    int k = 0;
    for (int b = 0; b < m->nb; b++) {
        float d = seg_dist(p, &m->bones[b]);
        if (d < CUTOFF) { w[b] = 1.0f / powf(d + WEIGHT_D, WEIGHT_P); k++; }
        else            { w[b] = 0.0f; }
    }
    if (k == 0) {   /* always bind to the nearest bone */
        int best = 0; float bd = 1e30f;
        for (int b = 0; b < m->nb; b++) {
            float d = seg_dist(p, &m->bones[b]);
            if (d < bd) { bd = d; best = b; }
        }
        w[best] = 1.0f; k = 1;
    }
    return k;
}

/* ordering key for one bone */
static float influence_score(float w, float total, float r, float R, float beta) {
    float wn = w / total;
    float d = r / R;
    if (d < 1e-6f) d = 1e-6f;
    return wn * powf(d, beta);
}

/* full and truncated skinning of one vertex; returns the displacement error */
static float skin_vertex(const Model* m, V3 p, float beta, Skin* s,
                         V3* full, V3* trunc) {
    float w[MAXB], total = 0.0f;
    memset(s, 0, sizeof(*s));
    s->ninfl = raw_weights(m, p, w);
    for (int b = 0; b < m->nb; b++) total += w[b];

    /* score every live influence, and record the weight factor */
    int order[MAXB]; int n = 0;
    for (int b = 0; b < m->nb; b++) {
        if (w[b] <= 0.0f) continue;
        float r = vlen(vsub(p, m->bones[b].p0));
        s->score[b] = influence_score(w[b], total, r, m->radius, beta);
        order[n++] = b;
    }
    /* insertion sort, descending score, bone index breaks ties */
    for (int i = 1; i < n; i++) {
        int key = order[i], j = i - 1;
        while (j >= 0 && (s->score[order[j]] < s->score[key] ||
                          (s->score[order[j]] == s->score[key] && order[j] > key))) {
            order[j+1] = order[j]; j--;
        }
        order[j+1] = key;
    }

    int keep = n < MAXK ? n : MAXK;
    float kept_raw = 0.0f;
    for (int i = 0; i < keep; i++) kept_raw += w[order[i]];
    s->nkept = keep;
    for (int i = 0; i < keep; i++) {
        int b = order[i];
        s->kept[i] = b;
        s->w[b] = w[b] / kept_raw;
    }
    s->loss = total > 0.0f ? (total - kept_raw) / total : 0.0f;

    /* deform with every influence, then with the kept ones only */
    V3 wf = v3(0,0,0), wt = v3(0,0,0);
    for (int b = 0; b < m->nb; b++) {
        if (w[b] <= 0.0f) continue;
        V3 d = rotz(p, m->bones[b].p0, m->bones[b].angle);
        wf = vadd(wf, vmul(d, w[b] / total));
        if (s->w[b] > 0.0f) wt = vadd(wt, vmul(d, s->w[b]));
    }
    if (full)  *full = wf;
    if (trunc) *trunc = wt;
    return vlen(vsub(wf, wt));
}

typedef struct {
    int over;
    float max_loss;
    float bump;         /* max displacement error, i.e. bump height */
    float rms;          /* rms displacement error, the solver objective */
} Stats;

static Stats measure(const Model* m, float beta) {
    Stats st = {0, 0.0f, 0.0f, 0.0f};
    double sq = 0.0;
    for (int i = 0; i < m->nv; i++) {
        Skin s;
        float e = skin_vertex(m, m->rest[i], beta, &s, NULL, NULL);
        if (s.ninfl > MAXK) st.over++;
        if (s.loss > st.max_loss) st.max_loss = s.loss;
        if (e > st.bump) st.bump = e;
        sq += (double)e * e;
    }
    st.rms = (float)sqrt(sq / (double)m->nv);
    return st;
}

/* back-solve beta from the measured bump height: the smallest beta that
 * attains the minimum maximum-truncation-displacement over the mesh.  Bump
 * is piecewise constant in beta, so the minimum is usually a plateau; report
 * both ends rather than an arbitrary interior point. */
static float solve_beta(const Model* m, float* out_bump, float* out_hi) {
    int n = (int)(BETA_MAX / BETA_STEP) + 1;
    static float bump[1024];
    if (n > (int)(sizeof(bump)/sizeof(bump[0]))) n = (int)(sizeof(bump)/sizeof(bump[0]));
    float best = 1e30f;
    for (int i = 0; i < n; i++) {
        bump[i] = measure(m, (float)i * BETA_STEP).bump;
        if (bump[i] < best) best = bump[i];
    }
    float lo = -1.0f, hi = -1.0f;
    for (int i = 0; i < n; i++) {
        if (bump[i] <= best + 1e-7f) {
            if (lo < 0.0f) lo = (float)i * BETA_STEP;
            hi = (float)i * BETA_STEP;
        }
    }
    if (out_bump) *out_bump = best;
    if (out_hi)   *out_hi = hi;
    return lo < 0.0f ? 0.0f : lo;
}

/* how the ordering key splits between its two factors, in log space:
 * |ln weight| and |beta ln(r/R)| over every live influence. */
static void factor_share(const Model* m, float beta, float* w_share,
                         float* d_share, float* ratio) {
    double lw = 0.0, ld = 0.0;
    for (int i = 0; i < m->nv; i++) {
        float w[MAXB], total = 0.0f;
        raw_weights(m, m->rest[i], w);
        for (int b = 0; b < m->nb; b++) total += w[b];
        for (int b = 0; b < m->nb; b++) {
            if (w[b] <= 0.0f) continue;
            double wn = (double)w[b] / (double)total;
            double r = (double)vlen(vsub(m->rest[i], m->bones[b].p0)) / (double)m->radius;
            if (r < 1e-6) r = 1e-6;
            lw += fabs(log(wn));
            ld += fabs((double)beta * log(r));
        }
    }
    double sum = lw + ld;
    if (sum <= 0.0) { *w_share = 1.0f; *d_share = 0.0f; *ratio = 0.0f; return; }
    *w_share = (float)(lw / sum);
    *d_share = (float)(ld / sum);
    *ratio = lw > 0.0 ? (float)(ld / lw) : 0.0f;   /* distance : weight */
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void self_test(void) {
    /* score ordering: a light bone on a long lever must beat a heavier bone
     * on a short lever once beta > 0 */
    float wA = 0.25f, wB = 0.20f, R = 2.0f;
    float sA0 = influence_score(wA, 1.0f, 1.0f, R, 0.0f);
    float sB0 = influence_score(wB, 1.0f, 2.0f, R, 0.0f);
    check("beta0_orders_by_weight", sA0 > sB0);
    float sA1 = influence_score(wA, 1.0f, 1.0f, R, 1.0f);
    float sB1 = influence_score(wB, 1.0f, 2.0f, R, 1.0f);
    check("beta1_orders_by_lever", sB1 > sA1);

    Model m;
    build_model(&m, "selftest", 6);
    float wr, dr, ratio;
    factor_share(&m, 0.0f, &wr, &dr, &ratio);
    check("beta0_no_distance_share", dr == 0.0f && wr == 1.0f);
    factor_share(&m, 1.0f, &wr, &dr, &ratio);
    check("beta1_split_sums_to_one", fabsf(wr + dr - 1.0f) < 1e-5f);
    check("beta1_has_distance_share", dr > 0.0f && dr < 1.0f);

    Stats st = measure(&m, 1.0f);
    check("over_limit_present", st.over > 0);
    check("loss_in_range", st.max_loss >= 0.0f && st.max_loss < 1.0f);
    check("bump_nonnegative", st.bump >= 0.0f);

    /* every vertex ends with at most four influences, renormalized */
    int bad_keep = 0, bad_norm = 0;
    for (int i = 0; i < m.nv; i++) {
        Skin s;
        skin_vertex(&m, m.rest[i], 1.0f, &s, NULL, NULL);
        if (s.nkept > MAXK) bad_keep++;
        float sum = 0.0f;
        for (int b = 0; b < m.nb; b++) sum += s.w[b];
        if (fabsf(sum - 1.0f) > 1e-5f) bad_norm++;
    }
    check("budget_respected", bad_keep == 0);
    check("weights_renormalized", bad_norm == 0);

    float bump;
    float beta = solve_beta(&m, &bump, NULL);
    check("beta_in_range", beta >= 0.0f && beta <= BETA_MAX);
    check("solved_not_worse_than_zero", bump <= measure(&m, 0.0f).bump + 1e-6f);
}

int main(int argc, char** argv) {
    int verbose = argc > 1 && strcmp(argv[1], "-v") == 0;
    self_test();
    printf("\n");

    /* Every group must actually truncate: with 2 axial bones plus the collar
     * ring, a model of density d has d bones in total, so d must exceed 4. */
    const int dens[3] = { 6, 8, 12 };
    Stats table[3];
    float betas[3], betas_hi[3], wshare[3], dshare[3];

    printf("%-8s %3s %5s %6s %9s %9s %9s %8s %9s %9s\n",
           "model", "nb", "over", "lossMax", "bumpWo", "bumpDist", "beta", "wShare", "dShare", "d:w");
    for (int g = 0; g < 3; g++) {
        Model m;
        char name[16];
        snprintf(name, sizeof(name), "density%d", dens[g]);
        build_model(&m, name, dens[g]);

        Stats wo = measure(&m, 0.0f);
        float bump, hi;
        float beta = solve_beta(&m, &bump, &hi);
        Stats st = measure(&m, beta);
        float wr, dr, rt;
        factor_share(&m, beta, &wr, &dr, &rt);

        table[g] = st; betas[g] = beta; betas_hi[g] = hi;
        wshare[g] = wr; dshare[g] = dr;

        if (st.over == 0) failures++;

        printf("%-8s %3d %5d %6.4f %9.6f %9.6f %9.4f %8.4f %9.4f %9.4f\n",
               m.name, m.nb, st.over, st.max_loss, wo.bump, st.bump,
               beta, wr, dr, rt);
    }

    /* the twelve required numbers */
    printf("\nmodel      overLimit  maxWeightLoss  bumpHeight  weightShare  distShare  betaLo  betaHi\n");
    for (int g = 0; g < 3; g++) {
        printf("density%-3d %8d  %11.6f  %10.6f  %11.4f  %9.4f  %6.3f  %6.3f\n",
               dens[g], table[g].over, table[g].max_loss, table[g].bump,
               wshare[g], dshare[g], betas[g], betas_hi[g]);
    }
    if (verbose) {
        printf("\nbump vs beta (bumpWo is beta=0):\n");
        printf("%-10s", "beta");
        for (int g = 0; g < 3; g++) printf("  density%-3d", dens[g]);
        printf("\n");
        for (int i = 0; i <= 30; i++) {
            float beta = (float)i * 0.2f;
            printf("%-10.2f", beta);
            for (int g = 0; g < 3; g++) {
                Model m;
                build_model(&m, "sweep", dens[g]);
                printf("  %9.6f", measure(&m, beta).bump);
            }
            printf("\n");
        }
    }

    int hi = 0;
    for (int g = 1; g < 3; g++)
        if (table[g].bump > table[hi].bump) hi = g;
    printf("\nhighest_bump=density%d\n", dens[hi]);
    if (hi == 2) {
        printf("verdict: highest-density bump is the largest -> distance term is too small.\n");
    } else {
        printf("verdict: highest-density bump is not the largest -> distance term is adequate.\n");
    }
    printf("  density%d: bump=%.6f  weightShare=%.4f  distanceShare=%.4f  betaRatio(d:w)=%.3f:1  betaRange=[%.2f,%.2f]\n",
           dens[hi], table[hi].bump, wshare[hi], dshare[hi], betas[hi], betas[hi], betas_hi[hi]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
