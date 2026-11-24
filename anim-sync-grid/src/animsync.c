/*
 * animsync: align an animation's start to a fixed time grid.
 *
 * Grid width G is derived jointly from the frame interval and the network
 * sync precision.  A width is adopted only when BOTH constraints hold:
 *
 *   (F) frame rate      G = k * Tf,  k a positive integer
 *                       -> every grid point falls on a frame boundary, so a
 *                          snapped start can always begin on a frame, and two
 *                          devices that land in the same cell show the same
 *                          frame; adjacent cells differ by exactly k frames.
 *   (N) network sync    G >= 2 * sigma
 *                       -> one cell spans the whole mutual clock-error band,
 *                          so devices that agree to within the precision
 *                          coalesce into the same cell instead of scattering
 *                          across many cells.
 *
 * Smallest width satisfying both:  k = ceil(2*sigma / Tf),  G = k * Tf.
 * A larger sigma pushes k up; a higher fps shrinks Tf and pushes k up as well.
 * That is the "width derived from both" rule.
 *
 * Playback model, per device i:
 *   clock estimate   r_i = t0 + eps_i,   eps_i in [-sigma, +sigma]
 *   grid snap        s_i = round(r_i / G) * G
 *   frame start      b_i = ceil(s_i / Tf) * Tf     (playback begins on a frame)
 *   frame at time t  f_i(t) = floor((t - b_i) / Tf)
 *   cross-device desync  max_t ( max_i f_i(t) - min_i f_i(t) )
 *
 * f_i is monotone in b_i, so for a pair the desync is ceil(|b_i - b_j| / Tf),
 * maximised at the eps extremes.  A tier is swept over every grid phase t0 in
 * [0, G), i.e. the worst case over repeated plays.
 *
 * Per frame rate we record (the "nine numbers"):
 *   start offset deviation : max_i |b_i - t0|                    [ms]
 *   grid width utilisation : 100 * 2*sigma / G                    [%]
 *                            how much of a cell the network band fills;
 *                            100% = cell sized purely by the network,
 *                            lower = the frame rate coarsened the cell.
 *   max cross-device desync: as defined above                    [frames]
 *
 * A network-only grid (G = 2*sigma, frame rate ignored) is evaluated too, to
 * show what the failure actually looks like.
 */
#include <stdio.h>
#include <math.h>

#define PHASES 20000

typedef struct {
    double fps, tf, sigma, band, g;
    int    k;            /* G / Tf */
    int    aligned;      /* G is an integer number of frames */
    int    network_only; /* G set by network precision alone */
} Tier;

static Tier derive_joint(double fps, double sigma) {
    Tier t;
    t.fps = fps;  t.tf = 1000.0 / fps;  t.sigma = sigma;  t.band = 2.0 * sigma;
    t.k = (int)ceil(t.band / t.tf - 1e-9);        /* both constraints */
    if (t.k < 1) t.k = 1;
    t.g = t.k * t.tf;
    t.aligned = 1;
    t.network_only = 0;
    return t;
}

static Tier derive_network_only(double fps, double sigma) {
    Tier t;
    t.fps = fps;  t.tf = 1000.0 / fps;  t.sigma = sigma;  t.band = 2.0 * sigma;
    t.g = t.band;                                  /* frame rate ignored */
    t.k = 0;
    t.aligned = (fabs(t.g / t.tf - round(t.g / t.tf)) < 1e-9);
    t.network_only = 1;
    return t;
}

typedef struct {
    double start_dev_ms;   /* max |b - t0|                */
    int    desync_fr;      /* max cross-device desync      */
    double desync_ms;      /* same, expressed in ms        */
} Measured;

static Measured measure(const Tier* t) {
    Measured m = { 0.0, 0, 0.0 };
    for (int p = 0; p < PHASES; p++) {
        double t0 = t->g * (double)p / (double)PHASES;
        double b_lo = 0.0, b_hi = 0.0;
        int first = 1;

        for (int e = 0; e < 2; e++) {
            double eps = (e == 0) ? -t->sigma : t->sigma;
            double r   = t0 + eps;
            double s   = round(r / t->g) * t->g;
            double b   = ceil(s / t->tf - 1e-9) * t->tf;

            double dev = fabs(b - t0);
            if (dev > m.start_dev_ms) m.start_dev_ms = dev;

            if (first) { b_lo = b_hi = b; first = 0; }
            else { if (b < b_lo) b_lo = b; if (b > b_hi) b_hi = b; }
        }
        double span = b_hi - b_lo;
        if (span > m.desync_ms) m.desync_ms = span;
        int fr = (int)ceil(span / t->tf - 1e-9);
        if (fr < 0) fr = 0;
        if (fr > m.desync_fr) m.desync_fr = fr;
    }
    return m;
}

/* ---------------- self-test scaffolding (repo style) ---------------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("check %-40s %s\n", name, ok ? "ok" : "FAIL");
    if (!ok) failures++;
}

int main(void) {
    const double sigma = 20.0;                 /* one-way clock sync precision */
    const double fps[3] = { 24.0, 30.0, 60.0 };

    Tier     tj[3], tn[3];
    Measured mj[3], mn[3];

    printf("net_sync_precision sigma = %.1f ms   (band 2*sigma = %.1f ms)\n",
           sigma, 2.0 * sigma);
    printf("joint rule: k = ceil(2*sigma/Tf), G = k*Tf  (frame-aligned AND >= band)\n\n");

    for (int i = 0; i < 3; i++) {
        tj[i] = derive_joint(fps[i], sigma);
        mj[i] = measure(&tj[i]);
        tn[i] = derive_network_only(fps[i], sigma);
        mn[i] = measure(&tn[i]);
    }

    /* ---- the nine numbers, joint derivation ---- */
    printf("%-6s %-9s %-6s %-9s %-8s %-15s %-9s %-11s\n",
           "fps", "Tf(ms)", "k", "G(ms)", "G/Tf", "start_dev(ms)",
           "util(%)", "desync(fr)");
    for (int i = 0; i < 3; i++) {
        printf("%-6.0f %-9.3f %-6d %-9.3f %-8.3f %-15.3f %-9.2f %-11d\n",
               tj[i].fps, tj[i].tf, tj[i].k, tj[i].g, tj[i].g / tj[i].tf,
               mj[i].start_dev_ms, 100.0 * tj[i].band / tj[i].g, mj[i].desync_fr);
    }

    printf("\nNINE_NUMBERS  fps : start_dev_ms  util_pct  desync_fr\n");
    for (int i = 0; i < 3; i++)
        printf("  %2.0f : %8.3f  %7.2f  %d\n",
               tj[i].fps, mj[i].start_dev_ms,
               100.0 * tj[i].band / tj[i].g, mj[i].desync_fr);

    int low_most = (mj[0].desync_fr > mj[1].desync_fr &&
                    mj[0].desync_fr > mj[2].desync_fr);
    int worst = 0;
    for (int i = 1; i < 3; i++)
        if (mj[i].desync_fr > mj[worst].desync_fr) worst = i;

    printf("\nlow_fps_has_most_desync_frames=%d   (lowest-fps tier has the FEWEST: %s)\n",
           low_most,
           (mj[0].desync_fr < mj[1].desync_fr && mj[0].desync_fr < mj[2].desync_fr)
               ? "yes" : "no");
    printf("most_desync_tier=%.0f fps\n", tj[worst].fps);

    /* ---- grid width vs frame interval, side by side ---- */
    printf("\nGRID_VS_FRAME\n");
    printf("%-6s %-12s %-12s %-13s %-8s\n",
           "fps", "frame(ms)", "grid(ms)", "grid/frames", "aligned");
    for (int i = 0; i < 3; i++)
        printf("%-6.0f %-12.3f %-12.3f %-13.2f %-8s\n",
               tj[i].fps, tj[i].tf, tj[i].g, tj[i].g / tj[i].tf,
               tj[i].aligned ? "yes" : "NO");

    /* ---- what network-only actually looks like ---- */
    printf("\nNETWORK_ONLY_GRID  G = 2*sigma, frame rate ignored\n");
    printf("%-6s %-9s %-10s %-11s %-12s\n",
           "fps", "G(ms)", "G/Tf", "aligned", "desync(fr)");
    for (int i = 0; i < 3; i++)
        printf("%-6.0f %-9.3f %-10.3f %-11s %-12d\n",
               tn[i].fps, tn[i].g, tn[i].g / tn[i].tf,
               tn[i].aligned ? "yes" : "NO", mn[i].desync_fr);
    printf("  (G is identical at every fps, and not frame-aligned; desync is"
           " no better at any tier -> this is the real signature of a"
           " network-only width, not 'low fps is worst')\n");

    /* ---- why the stated diagnostic cannot fire: desync_fr = span/Tf ---- */
    printf("\nSIGMA_SWEEP  desync_fr per fps   (frame-count desync scales as 1/Tf)\n");
    printf("%-10s %-8s %-8s %-8s\n", "sigma(ms)", "24fps", "30fps", "60fps");
    const double sw[6] = { 4.0, 8.0, 12.0, 20.0, 30.0, 40.0 };
    int ever_low_strict = 0;
    for (int s = 0; s < 6; s++) {
        int df[3];
        for (int i = 0; i < 3; i++) {
            Tier tt = derive_joint(fps[i], sw[s]);
            df[i] = measure(&tt).desync_fr;
        }
        if (df[0] > df[1] && df[0] > df[2]) ever_low_strict = 1;
        printf("%-10.1f %-8d %-8d %-8d\n", sw[s], df[0], df[1], df[2]);
    }

    printf("\n");
    for (int i = 0; i < 3; i++) {
        char name[64];
        snprintf(name, sizeof name, "tier_%.0f_covers_band", tj[i].fps);
        check(name, tj[i].g >= tj[i].band - 1e-9);
        snprintf(name, sizeof name, "tier_%.0f_frame_aligned", tj[i].fps);
        check(name, fabs(tj[i].g / tj[i].tf - tj[i].k) < 1e-9);
        snprintf(name, sizeof name, "tier_%.0f_desync_eq_k", tj[i].fps);
        check(name, mj[i].desync_fr == tj[i].k);
    }
    check("network_only_same_g_all_fps",
          fabs(tn[0].g - tn[1].g) < 1e-9 && fabs(tn[1].g - tn[2].g) < 1e-9);
    check("network_only_not_aligned_at_low_fps", !tn[0].aligned);
    check("joint_beats_network_only_at_low_fps",
          mj[0].desync_fr <= mn[0].desync_fr);
    check("desync_fr_strictly_rises_with_fps_at_sigma20",
          mj[0].desync_fr < mj[1].desync_fr && mj[1].desync_fr < mj[2].desync_fr);
    check("low_fps_never_strictly_worst_over_sweep", !ever_low_strict);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
