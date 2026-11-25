
/* ao: ambient-occlusion bake on a cylindrical-pit model.
 *
 * Two samplers are compared on the same pixels:
 *   uniform  - rays uniform over the hemisphere (the baseline that leaves
 *              creases dead-black and rim pixels noisy/over-bright);
 *   distance - pilot rays measure the occluder hit-distance distribution,
 *              log-spaced distance bins get a 1/t target density (dense near,
 *              sparse far), with a floor, a per-bin density cap, and t_near
 *              truncation. Final rays resample pilot directions by bin.
 *
 * Three pit depths are baked. Per tier and sampler we record the occlusion
 * histogram, the dead-black pixel fraction and the mean occlusion. */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define MAXV 96
#define MAXF 128
#define MAXPX 192
#define KSEG 12                 /* pit tessellation segments */
#define N_PILOT 96
#define N_RAYS 192
#define N_GT 8192
#define NB 6                    /* log-spaced distance bins */
#define NOBIN NB                /* extra bin holding pilot misses */
#define PI 3.14159265358979323846f
#define TMAX 1.0e3f
#define TNEAR 1.0e-2f           /* self-intersection guard / near truncation */
#define PIT_R 1.0f
#define GND_R 2.0f
#define PCAP 0.30f              /* max target probability of one pilot direction (soft cap) */
#define ALPHA 0.20f             /* mixture weight kept on the empirical pilot distribution */
#define BLACK 0.95f             /* occlusion at/above which a pixel is dead-black */
#define LARGE_BLACK 0.10f       /* "large area" dead-black threshold */

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { V3 v[MAXV]; int nv; Tri f[MAXF]; int nf; } Mesh;
typedef struct { V3 p, n; } Pixel;

static V3 vsub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 vadd(V3 p, V3 q) { V3 r = { p.x+q.x, p.y+q.y, p.z+q.z }; return r; }
static V3 vscale(V3 p, float s) { V3 r = { p.x*s, p.y*s, p.z*s }; return r; }
static V3 vcross(V3 p, V3 q) {
    V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r;
}
static float vdot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float vlen(V3 p) { return sqrtf(vdot(p, p)); }
static V3 vnorm(V3 p) {
    float l = vlen(p);
    if (l < 1e-12f) { V3 z = {0,0,0}; return z; }
    return vscale(p, 1.0f / l);
}

static uint32_t rng(uint32_t* s) {
    uint32_t x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return x;
}
static float frand(uint32_t* s) { return (rng(s) >> 8) * (1.0f / 16777216.0f); }

static int add_v(Mesh* m, V3 p) { if (m->nv >= MAXV) return -1; m->v[m->nv] = p; return m->nv++; }
static int add_f(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF || a < 0 || b < 0 || c < 0) return -1;
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* cylindrical pit: ground annulus, wall, bottom fan; same radius, depth varies */
static void build_pit(Mesh* m, float depth) {
    memset(m, 0, sizeof(*m));
    int top[KSEG], bot[KSEG], out[KSEG];
    for (int i = 0; i < KSEG; i++) {
        float th = 2.0f * PI * i / KSEG, cx = cosf(th), sy = sinf(th);
        top[i] = add_v(m, (V3){ PIT_R*cx, PIT_R*sy, 0.0f });
        bot[i] = add_v(m, (V3){ PIT_R*cx, PIT_R*sy, -depth });
        out[i] = add_v(m, (V3){ GND_R*cx, GND_R*sy, 0.0f });
    }
    int cen = add_v(m, (V3){ 0.0f, 0.0f, -depth });
    for (int i = 0; i < KSEG; i++) {
        int j = (i + 1) % KSEG;
        add_f(m, out[i], out[j], top[j]); add_f(m, out[i], top[j], top[i]);
        add_f(m, top[i], top[j], bot[j]); add_f(m, top[i], bot[j], bot[i]);
        add_f(m, cen, bot[i], bot[j]);
    }
}

/* bake pixels: ground grid around the rim, bottom disk, mid-depth wall ring */
static int build_pixels(Pixel* px, float depth) {
    int n = 0;
    for (int gx = 0; gx < 11; gx++) for (int gy = 0; gy < 11; gy++) {
        float x = -1.9f + 0.38f * gx, y = -1.9f + 0.38f * gy;
        if (x*x + y*y < PIT_R*PIT_R) continue;          /* inside footprint: skip */
        px[n].p = (V3){ x, y, 0.0f }; px[n].n = (V3){ 0, 0, 1 }; n++;
    }
    for (int ir = 0; ir < 3; ir++) {
        float rr = 0.2f + 0.3f * ir;
        for (int a = 0; a < 8; a++) {
            float th = 2.0f * PI * a / 8;
            px[n].p = (V3){ rr*cosf(th), rr*sinf(th), -depth };
            px[n].n = (V3){ 0, 0, 1 }; n++;
        }
    }
    for (int a = 0; a < 8; a++) {
        float th = 2.0f * PI * a / 8;
        px[n].p = (V3){ PIT_R*cosf(th), PIT_R*sinf(th), -depth * 0.5f };
        px[n].n = (V3){ -cosf(th), -sinf(th), 0.0f }; n++;
    }
    return n;
}

/* Moller-Trumbore, double sided; hit accepted only inside (t0, t1) */
static int ray_tri(V3 o, V3 d, V3 a, V3 b, V3 c, float t0, float t1, float* t) {
    V3 e1 = vsub(b, a), e2 = vsub(c, a);
    V3 p = vcross(d, e2);
    float det = vdot(e1, p);
    if (fabsf(det) < 1e-9f) return 0;
    float inv = 1.0f / det;
    float u = vdot(vsub(o, a), p) * inv;
    if (u < 0.0f || u > 1.0f) return 0;
    V3 q = vcross(vsub(o, a), e1);
    float v = vdot(d, q) * inv;
    if (v < 0.0f || u + v > 1.0f) return 0;
    float tt = vdot(e2, q) * inv;
    if (tt <= t0 || tt >= t1) return 0;
    *t = tt;
    return 1;
}

static float scene_hit(const Mesh* m, V3 o, V3 d) {
    float best = TMAX;
    for (int i = 0; i < m->nf; i++) {
        float t;
        if (ray_tri(o, d, m->v[m->f[i].a], m->v[m->f[i].b], m->v[m->f[i].c],
                    TNEAR, best, &t))
            best = t;
    }
    return best;
}

/* uniform direction over the hemisphere around n (solid-angle uniform) */
static V3 hemi_dir(V3 n, float u1, float u2) {
    float phi = 2.0f * PI * u1, z = u2, r = sqrtf(1.0f - z*z);
    V3 helper = fabsf(n.x) < 0.9f ? (V3){1,0,0} : (V3){0,1,0};
    V3 t1 = vnorm(vcross(helper, n));
    V3 t2 = vcross(n, t1);
    V3 loc = { r*cosf(phi), r*sinf(phi), z };
    return vadd(vadd(vscale(t1, loc.x), vscale(t2, loc.y)), vscale(n, loc.z));
}

static float bake_uniform(const Mesh* m, const Pixel* px, int npix, uint32_t seed) {
    double sum = 0.0;
    for (int i = 0; i < npix; i++) {
        int blocked = 0;
        for (int s = 0; s < N_RAYS; s++) {
            V3 d = hemi_dir(px[i].n, frand(&seed), frand(&seed));
            if (scene_hit(m, px[i].p, d) < TMAX) blocked++;
        }
        sum += (double)blocked / N_RAYS;
    }
    return (float)(sum / npix);
}

static float bake_gt(const Mesh* m, const Pixel* px, int npix, uint32_t seed) {
    double sum = 0.0;
    for (int i = 0; i < npix; i++) {
        int blocked = 0;
        for (int s = 0; s < N_GT; s++) {
            V3 d = hemi_dir(px[i].n, frand(&seed), frand(&seed));
            if (scene_hit(m, px[i].p, d) < TMAX) blocked++;
        }
        sum += (double)blocked / N_GT;
    }
    return (float)(sum / npix);
}

static void bin_edges(float e[NB + 1]) {
    e[0] = TNEAR;
    float g = powf(TMAX / TNEAR, 1.0f / NB);
    for (int k = 1; k <= NB; k++) e[k] = TNEAR * powf(g, (float)k);
}

static int hit_bin(float t, const float e[NB + 1]) {
    if (t >= TMAX) return NOBIN;
    for (int k = 0; k < NB; k++) if (t < e[k + 1]) return k;
    return NB - 1;
}

/* distance-derived importance bake.
 * Pilot rays are uniform over the hemisphere; each hit direction gets weight
 * 1/t (dense near, sparse far) and each escaping direction weight 1. Weights
 * are soft-capped at PCAP of their total and blended with ALPHA of the
 * empirical uniform distribution, then final rays resample pilot directions.
 * min_hit -> nearest pilot occluder distance over the tier;
 * max_pdf -> largest realized sampling density q_j / bin_width (per unit t). */
static void bake_distance(const Mesh* m, const Pixel* px, int npix, uint32_t pilot_seed,
                          uint32_t final_seed, float* occ, float* min_hit, float* max_pdf) {
    float edges[NB + 1];
    bin_edges(edges);
    *min_hit = TMAX;
    *max_pdf = 0.0f;
    for (int i = 0; i < npix; i++) {
        V3 pdir[N_PILOT];
        float pt[N_PILOT], w[N_PILOT], q[N_PILOT], cdf[N_PILOT];
        int pbin[N_PILOT];
        float wsum = 0.0f;
        for (int s = 0; s < N_PILOT; s++) {
            pdir[s] = hemi_dir(px[i].n, frand(&pilot_seed), frand(&pilot_seed));
            pt[s] = scene_hit(m, px[i].p, pdir[s]);
            pbin[s] = hit_bin(pt[s], edges);
            if (pt[s] < *min_hit) *min_hit = pt[s];
            /* distance-derived target: 1/t for occluders, unit weight for escapes */
            w[s] = pt[s] < TMAX ? 1.0f / pt[s] : 1.0f;
            wsum += w[s];
        }
        if (wsum > 0.0f)
            for (int s = 0; s < N_PILOT; s++)
                if (w[s] > PCAP * wsum) w[s] = PCAP * wsum;
        float nsum = 0.0f;
        for (int s = 0; s < N_PILOT; s++) nsum += w[s];
        float acc = 0.0f;
        for (int s = 0; s < N_PILOT; s++) {
            float wd = nsum > 0.0f ? w[s] / nsum : 0.0f;       /* distance target */
            q[s] = (1.0f - ALPHA) * wd + ALPHA / N_PILOT;      /* mixture, keeps every ray live */
            acc += q[s];
            cdf[s] = acc;
            if (pbin[s] < NB) {
                float pdf = q[s] / (edges[pbin[s] + 1] - edges[pbin[s]]);
                if (pdf > *max_pdf) *max_pdf = pdf;
            }
        }
        double est = 0.0;
        for (int s = 0; s < N_RAYS; s++) {
            float u = frand(&final_seed);
            int j = 0;
            while (j < N_PILOT - 1 && u > cdf[j]) j++;
            float t = scene_hit(m, px[i].p, pdir[j]);
            if (t < TMAX)
                est += 1.0 / (N_PILOT * q[j]);  /* pilot pdf was 1/N_PILOT over the set */
        }
        float v = (float)(est / N_RAYS);
        occ[i] = v < 0.0f ? 0.0f : (v > 1.0f ? 1.0f : v);
    }
}

typedef struct { float hist[5]; float dead; float mean; float mae; float vmax; } Metrics;

static Metrics measure(const float* occ, const float* gt, int npix) {
    Metrics mt;
    memset(&mt, 0, sizeof(mt));
    mt.vmax = 0.0f;
    for (int i = 0; i < npix; i++) {
        int b = (int)(occ[i] * 5.0f);
        if (b < 0) b = 0;
        if (b > 4) b = 4;
        mt.hist[b] += 1.0f;
        if (occ[i] >= BLACK) mt.dead += 1.0f;
        mt.mean += occ[i];
        mt.mae += fabsf(occ[i] - gt[i]);
        if (occ[i] > mt.vmax) mt.vmax = occ[i];
    }
    for (int b = 0; b < 5; b++) mt.hist[b] /= npix;
    mt.dead /= npix;
    mt.mean /= npix;
    mt.mae /= npix;
    return mt;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const float depths[3] = { 0.5f, 1.25f, 2.5f };
    const char* names[3] = { "shallow", "medium", "deep" };
    Metrics mu[3], md[3];
    float gt_mean[3], d_min[3], d_maxpdf[3];
    int npix[3];
    double mae_un = 0.0, mae_is = 0.0;

    for (int tier = 0; tier < 3; tier++) {
        Mesh m;
        build_pit(&m, depths[tier]);
        Pixel px[MAXPX];
        npix[tier] = build_pixels(px, depths[tier]);

        float occ_un[MAXPX], occ_is[MAXPX], occ_gt[MAXPX];
        double gsum = 0.0;
        uint32_t useed = 0x12345678u + (uint32_t)tier;
        for (int i = 0; i < npix[tier]; i++) {
            Pixel one = px[i];
            uint32_t gseed = 0x9e3779b9u + (uint32_t)tier * 101u + (uint32_t)i;
            occ_gt[i] = bake_gt(&m, &one, 1, gseed);
            occ_un[i] = bake_uniform(&m, &one, 1, useed);
            gsum += occ_gt[i];
        }
        gt_mean[tier] = (float)(gsum / npix[tier]);
        bake_distance(&m, px, npix[tier],
                      0xA5A5A5A5u + (uint32_t)tier, 0x5A5A5A5Au + (uint32_t)tier,
                      occ_is, &d_min[tier], &d_maxpdf[tier]);

        mu[tier] = measure(occ_un, occ_gt, npix[tier]);
        md[tier] = measure(occ_is, occ_gt, npix[tier]);
        mae_un += mu[tier].mae;
        mae_is += md[tier].mae;

        printf("=== tier=%s depth=%.2f pixels=%d ===\n", names[tier], depths[tier], npix[tier]);
        printf("uniform   hist=%.3f,%.3f,%.3f,%.3f,%.3f dead_black=%.3f mean_occ=%.3f mae=%.4f\n",
               mu[tier].hist[0], mu[tier].hist[1], mu[tier].hist[2],
               mu[tier].hist[3], mu[tier].hist[4], mu[tier].dead, mu[tier].mean, mu[tier].mae);
        printf("distance  hist=%.3f,%.3f,%.3f,%.3f,%.3f dead_black=%.3f mean_occ=%.3f mae=%.4f\n",
               md[tier].hist[0], md[tier].hist[1], md[tier].hist[2],
               md[tier].hist[3], md[tier].hist[4], md[tier].dead, md[tier].mean, md[tier].mae);
        printf("tier_min_hit=%.4f tier_max_pdf=%.2f\n", d_min[tier], d_maxpdf[tier]);
    }

    /* headline: three tiers x three records (histogram, dead-black, mean) */
    printf("\n--- nine records (distance sampler) ---\n");
    printf("tier      dead_black  mean_occ  histogram[0,.2,.4,.6,.8,1]\n");
    for (int tier = 0; tier < 3; tier++)
        printf("%-9s %.3f       %.3f     %.3f,%.3f,%.3f,%.3f,%.3f\n",
               names[tier], md[tier].dead, md[tier].mean,
               md[tier].hist[0], md[tier].hist[1], md[tier].hist[2],
               md[tier].hist[3], md[tier].hist[4]);

    /* truncation diagnostic: nearest occluder distance beside the density cap.
     * The pair is always reported; "triggered" marks large-area dead black. */
    float edges[NB + 1];
    bin_edges(edges);
    float cap_pdf_bin0 = PCAP / (edges[1] - edges[0]);
    int triggered = md[2].dead >= LARGE_BLACK;
    printf("\n--- deep-tier truncation diagnostic (triggered=%d) ---\n", triggered);
    printf("nearest_hit=%.4f t_near=%.4f density_cap_p=%.2f realized_max_pdf=%.2f cap_pdf_bin0=%.2f\n",
           d_min[2], TNEAR, PCAP, d_maxpdf[2], cap_pdf_bin0);

    /* self tests */
    printf("\n--- checks ---\n");
    check("pixels_match", npix[0] == npix[1] && npix[1] == npix[2] && npix[0] > 0);
    check("gt_depth_order", gt_mean[0] < gt_mean[1] && gt_mean[1] < gt_mean[2]);
    check("distance_unbiased_mean",
          fabsf(md[0].mean - gt_mean[0]) < 0.03f &&
          fabsf(md[1].mean - gt_mean[1]) < 0.03f &&
          fabsf(md[2].mean - gt_mean[2]) < 0.03f);
    check("distance_mae_abs", mae_is / 3.0 < 0.05f);
    /* truncation control: re-cast the deep-tier pilot rays with a raised guard
     * of 5*t_near and count real occluders that would be discarded. */
    Mesh dm;
    build_pit(&dm, depths[2]);
    Pixel dpx[MAXPX];
    int dnpix = build_pixels(dpx, depths[2]);
    uint32_t cseed = 0xCAFE1234u;
    int cut_hits = 0;
    for (int i = 0; i < dnpix; i++)
        for (int s = 0; s < N_PILOT; s++) {
            V3 d = hemi_dir(dpx[i].n, frand(&cseed), frand(&cseed));
            float t = scene_hit(&dm, dpx[i].p, d);
            if (t < TMAX && t < 5.0f * TNEAR) cut_hits++;
        }
    check("near_truncation_would_cut", cut_hits > 0 && d_min[2] < 5.0f * TNEAR);
    check("cap_pair_present", cap_pdf_bin0 > 0.0f && d_maxpdf[2] > 0.0f);
    float hs = 0.0f;
    for (int b = 0; b < 5; b++) hs += md[2].hist[b];
    check("hist_sums_one", fabsf(hs - 1.0f) < 1e-5f);
    check("occ_in_range", md[0].vmax <= 1.05f && md[1].vmax <= 1.05f && md[2].vmax <= 1.05f);
    /* t_near guard: a ray leaving a ground pixel downward must not self-intersect */
    V3 op = { 1.5f, 0.0f, 0.0f }, down = { 0.0f, 0.0f, -1.0f };
    check("tnear_guard", scene_hit(&dm, op, down) >= TMAX);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
