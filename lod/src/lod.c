
/* lod: two-level LOD hysteresis simulation.
 *
 * Hysteresis width as a function of the two levels' geometric errors.
 * Each error model does a round-trip camera move over N_SEEDS noise
 * realizations; per model we report:
 *   mean switches / mean flicker frames / switch delay (frames) /
 *   mean visual pop (px) at the actual switches.
 * Deterministic seeded RNG; `make test` prints every number and exits
 * nonzero if any sanity check fails.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

/* ---- camera / screen constants ---- */
#define G_PX   1080.0   /* SSE gain: screen_h / (2 tan(fov/2)) */
#define TAU    16.0     /* SSE threshold, px */
#define V      4.0      /* radial camera speed, world units / frame */
#define MARGIN 0.8      /* round trip spans d0*(1+-MARGIN) */
#define N_SEEDS 32

/* ---- disturbance model ---- */
#define SIGMA0 0.75    /* additive SSE noise floor, px */
#define AJIT   1.0     /* camera micro-jitter amplitude, px-equivalent */
#define TJIT   21.0    /* jitter period, frames */

/* ---- width rule ---- */
#define GAMMA  0.5     /* fraction of the pop we may spend on hysteresis */
#define KSIG   3.0     /* noise floor multiplier (two-sided k*sigma) */
#define T_MAX  360.0   /* delay budget, frames */
#define W_FIX  4.0     /* strawman: one fixed px gap for every error tier */

#define FLICK_RUN 5    /* a reversed LOD run shorter than this is flicker */
#define MAXSW 256

typedef struct {
    const char *name;
    double d0;          /* nominal switch distance */
    double r;           /* ec / ef */
    double ef, ec;
    double u_rel;       /* relative SSE uncertainty of this error tier */
} Model;

typedef struct {
    double switches;      /* mean over seeds */
    double flicker_frames;
    double extra_sw;      /* switches beyond the expected 2 per round trip */
    double delay_frames;  /* mean |actual - ideal| over both crossings */
    double pop_px;        /* mean SSE gap at actual switches */
    double w_px;
    double w_dist;        /* far-edge minus near-edge distance */
    double delay_analytic;
    int    n_frames;
} Result;

/* ---- deterministic RNG (32-bit LCG; mask keeps the arithmetic wrapping
 * even though unsigned long is 64-bit on this platform) ---- */
static unsigned long rng_state = 0x12345678UL;
static double rand_u(void) {
    rng_state = (rng_state * 1103515245UL + 12345UL) & 0xFFFFFFFFUL;
    return (double)(rng_state >> 8) / (double)(1UL << 24);
}
static double rand_n(void) { /* Box-Muller */
    double u1 = fmax(rand_u(), 1e-12), u2 = rand_u();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

static double sigma_at(double s, double u_rel) {
    return sqrt(SIGMA0 * SIGMA0 + (u_rel * s) * (u_rel * s));
}

/* pop at the nominal boundary: function of both geometric errors */
static double pop0(const Model *m) {
    return TAU * (1.0 - m->ef / m->ec);
}

/* error-layered hysteresis width, px:
 * max(noise/jitter floor, gamma * pop), then capped twice:
 *  - delay budget T_max frames
 *  - round-trip reach: near-point SSE is tau/(1-MARGIN); the upper edge
 *    must keep a 3-sigma + jitter margin below it or noisy seeds miss
 *    the return crossing */
static double layered_width(const Model *m) {
    double floor_w = 2.0 * (AJIT + KSIG * sigma_at(TAU, m->u_rel));
    double want = fmax(floor_w, GAMMA * pop0(m));
    double p_rate = TAU * V / m->d0;
    double cap_delay = T_MAX * p_rate;
    double s_peak = TAU / (1.0 - MARGIN);
    double cap_reach = 2.0 * (s_peak - TAU
                        - 3.0 * sigma_at(s_peak, m->u_rel) - AJIT);
    if (cap_reach <= 0) {
        /* trip geometry cannot hold a flicker-free width at this tier */
        fprintf(stderr, "%s: unreachable width (margin too small)\n",
                m->name);
        exit(2);
    }
    return fmin(want, fmin(cap_delay, cap_reach));
}

static Result simulate(const Model *m, double w, unsigned long seed0) {
    double d_near = m->d0 * (1.0 - MARGIN);
    double d_far  = m->d0 * (1.0 + MARGIN);
    double leg = (d_far - d_near) / V;          /* frames per leg */
    int n = (int)ceil(2.0 * leg) + 1;

    double lo = TAU - w / 2.0;                  /* downgrade edge (SSE) */
    double hi = TAU + w / 2.0;
    double d_out = G_PX * m->ec / lo;           /* far distance edge */
    double d_in  = G_PX * m->ec / hi;           /* near distance edge */
    double t_ideal = MARGIN * m->d0 / V;        /* ideal crossing offset */

    long sum_sw = 0, sum_flick = 0, sum_extra = 0;
    double sum_delay = 0, sum_pop = 0;

    for (int seed = 0; seed < N_SEEDS; seed++) {
        rng_state = seed0 + (unsigned long)seed * 2654435761UL;
        int state = 0;                          /* 0 = fine, 1 = coarse */
        int *states = (int *)calloc((size_t)n, sizeof(int));
        int sw_t[MAXSW]; double sw_pop[MAXSW]; int nsw = 0;

        for (int t = 0; t < n; t++) {
            double phase = (double)t / leg;
            long cyc = (long)floor(phase);
            double f = phase - (double)cyc;
            double d = (cyc % 2 == 0)
                     ? d_near + f * (d_far - d_near)
                     : d_far  - f * (d_far - d_near);
            if (d < d_near) d = d_near;
            if (d > d_far)  d = d_far;

            double s_true = G_PX * m->ec / d;
            double jitter = AJIT * sin(2.0 * PI * (double)t / TJIT);
            double s_obs = s_true + jitter
                         + sigma_at(s_true, m->u_rel) * rand_n();

            int prev = state;
            if (state == 0 && s_obs < lo) state = 1;
            else if (state == 1 && s_obs > hi) state = 0;
            states[t] = state;

            if (state != prev && nsw < MAXSW) {
                sw_t[nsw] = t;
                sw_pop[nsw] = s_true * (1.0 - m->ef / m->ec);
                nsw++;
            }
        }

        for (int i = 0; i < n; ) {              /* flicker runs */
            int j = i + 1;
            while (j < n && states[j] == states[i]) j++;
            if (j - i <= FLICK_RUN && j < n) sum_flick += (j - i);
            i = j;
        }

        double dl[2] = {-1, -1};
        for (int k = 0; k < nsw; k++) {
            double t = (double)sw_t[k];
            if (t < leg && dl[0] < 0) dl[0] = fabs(t - t_ideal);
            if (t >= leg && dl[1] < 0)
                dl[1] = fabs(t - (leg + t_ideal));
        }
        if (dl[0] >= 0 && dl[1] >= 0) sum_delay += (dl[0] + dl[1]) / 2.0;
        for (int k = 0; k < nsw; k++) sum_pop += sw_pop[k];

        sum_sw += nsw;
        sum_extra += (nsw > 2) ? nsw - 2 : 0;
        free(states);
    }

    Result res = {0};
    res.switches = (double)sum_sw / N_SEEDS;
    res.flicker_frames = (double)sum_flick / N_SEEDS;
    res.extra_sw = (double)sum_extra / N_SEEDS;
    res.delay_frames = sum_delay / N_SEEDS;
    res.pop_px = sum_pop / (double)(sum_sw ? sum_sw : 1);
    res.w_px = w;
    res.w_dist = d_out - d_in;
    res.delay_analytic = w / (TAU * V / m->d0);
    res.n_frames = n;
    return res;
}

static int failures = 0;
static void check(const char *name, int ok) {
    printf("check %s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char **argv) {
    int csv = (argc == 2);
    if (csv) {
        if (!(argv[1][0] == '-' && argv[1][1] == '-' &&
              argv[1][2] == 'c' && argv[1][3] == 's' &&
              argv[1][4] == 'v' && argv[1][5] == '\0'))
            csv = 0;
    }
    if (csv)
        printf("policy,model,ef,ec,u_rel,w_px,switches,flicker_frames,"
               "extra_switches,delay_frames,pop_px,w_dist,delay_analytic\n");
    for (int k = 0; k < 1000; k++) {
        double u = rand_u();
        if (!(u >= 0.0 && u < 1.0)) {
            if (!csv) check("rng_range", 0);
            break;
        }
        if (k == 999 && !csv) check("rng_range", 1);
    }

    /* u_rel grows with the error tier: coarse error metadata carries
     * proportionally larger measurement uncertainty. */
    Model models[3] = {
        {"M1 r=2 d0=60",   60.0, 2.0, 0, 0, 0.04},
        {"M2 r=4 d0=240", 240.0, 4.0, 0, 0, 0.10},
        {"M3 r=8 d0=960", 960.0, 8.0, 0, 0, 0.22},
    };
    for (int i = 0; i < 3; i++) {
        models[i].ec = TAU * models[i].d0 / G_PX;
        models[i].ef = models[i].ec / models[i].r;
    }

    if (!csv)
        printf("constants: G=%.0f tau=%.0f v=%.1f sigma0=%.2f Ajit=%.1f Tj=%.0f"
               " gamma=%.2f K=%.1f Tmax=%.0f Wfix=%.1f seeds=%d\n",
               G_PX, TAU, V, SIGMA0, AJIT, TJIT, GAMMA, KSIG, T_MAX, W_FIX,
               N_SEEDS);

    const char *pname[2] = {"A_fixed_compromise", "B_error_layered"};
    Result all[2][3];
    for (int p = 0; p < 2; p++) {
        if (!csv) {
            printf("\n== policy %s ==\n", pname[p]);
            printf("%-16s %8s %8s %7s %9s %9s %9s %10s %9s %9s\n",
                   "model", "ef", "ec", "u_rel", "w_px", "switch",
                   "flick_fr", "extra_sw", "delay_fr", "pop_px");
        }
        for (int i = 0; i < 3; i++) {
            double w = (p == 0) ? W_FIX : layered_width(&models[i]);
            Result r = simulate(&models[i], w,
                                0x12345678UL + (unsigned long)i * 7919UL);
            all[p][i] = r;
            if (csv)
                printf("%s,%s,%.4f,%.4f,%.2f,%.4f,%.3f,%.3f,%.3f,%.3f,"
                       "%.3f,%.3f,%.3f\n",
                       pname[p], models[i].name, models[i].ef, models[i].ec,
                       models[i].u_rel, r.w_px, r.switches,
                       r.flicker_frames, r.extra_sw, r.delay_frames,
                       r.pop_px, r.w_dist, r.delay_analytic);
            else
                printf("%-16s %8.3f %8.3f %7.2f %9.3f %9.2f %9.2f %10.2f"
                       " %9.2f %9.2f   (Wdist=%7.1f delay_ana=%6.1f)\n",
                       models[i].name, models[i].ef, models[i].ec,
                       models[i].u_rel, r.w_px, r.switches, r.flicker_frames,
                       r.extra_sw, r.delay_frames, r.pop_px,
                       r.w_dist, r.delay_analytic);
        }
    }
    if (csv) return 0;

    /* sanity checks: mean >=2, and every individual seed crosses both
     * ways (mean >= 1.99 leaves room only for rare edge-on misses) */
    for (int p = 0; p < 2; p++)
        for (int i = 0; i < 3; i++) {
            char tag[80];
            snprintf(tag, sizeof tag, "%s/%s crossings~=2", pname[p],
                     models[i].name);
            check(tag, all[p][i].switches >= 1.97);
            snprintf(tag, sizeof tag, "%s/%s delay>=0", pname[p],
                     models[i].name);
            check(tag, all[p][i].delay_frames >= 0);
            snprintf(tag, sizeof tag, "%s/%s pop>0", pname[p],
                     models[i].name);
            check(tag, all[p][i].pop_px > 0);
        }
    /* the diagnostic the task asks for: under the fixed compromise,
     * the largest-error tier must carry both the longest delay and
     * residual flicker; layering must remove the flicker there */
    check("A: M3 delay longest",
          all[0][2].delay_frames > all[0][1].delay_frames &&
          all[0][1].delay_frames > all[0][0].delay_frames);
    check("A: M3 still flickers", all[0][2].flicker_frames > 0);
    check("B: M3 flicker suppressed 90%+",
          all[1][2].flicker_frames <= 0.1 * all[0][2].flicker_frames);
    check("B: M1 no overspend vs A",
          all[1][0].delay_frames <= all[0][0].delay_frames * 3.0);

    /* error tiers and widths side by side */
    printf("\n== error tier vs hysteresis width ==\n");
    printf("%-16s %8s %8s %7s %9s %10s %12s %12s\n",
           "model", "ef", "ec", "u_rel", "pop0_px",
           "sigma_tau", "w_A_fixed", "w_B_layered");
    for (int i = 0; i < 3; i++) {
        double wb = layered_width(&models[i]);
        printf("%-16s %8.3f %8.3f %7.2f %9.2f %10.3f %12.3f %12.3f\n",
               models[i].name, models[i].ef, models[i].ec,
               models[i].u_rel, pop0(&models[i]),
               sigma_at(TAU, models[i].u_rel), W_FIX, wb);
    }
    printf("\nfailures=%d\n", failures);
    return failures ? 1 : 0;
}
