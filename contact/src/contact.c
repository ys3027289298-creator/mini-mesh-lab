/* contact: temporal smoothing of contact normals, keyed per contact surface.
 *
 * Scenario (2D): a disc slides down a 15-degree slope. The collider jitters
 * between frames, so the raw contact normal changes every frame (+-3 deg).
 * Without smoothing the body gets tangential kicks, slides along the surface
 * and takes a long time to be declared at rest.
 *
 * Fix: per-contact-surface temporal smoothing of the normal. The smoothing
 * coefficient is derived from the contact's own duration (age in frames):
 *
 *     alpha = 1 / (1 + age)
 *
 * A new contact has age 0 -> alpha 1 -> the raw normal is used as-is (new
 * contact points are NOT smoothed). The longer the contact lasts, the
 * smaller alpha gets and the heavier the smoothing.
 *
 * Experiment: three contact durations (short/medium/long, via initial
 * tangential speed). Per run we record:
 *   - normal angle variation: sum of frame-to-frame angle changes of the
 *     normal actually used for the response, in degrees
 *   - slip distance: path length travelled along the surface while in
 *     contact, in meters
 *   - frames until rest is declared (|v| < 0.02 for 5 consecutive frames)
 *
 * Modes:
 *   raw    - no smoothing (baseline)
 *   age    - alpha from the contact's own duration (the fix)
 *   global - alpha from the global frame count (the bug signature: the
 *            coefficient no longer reflects this contact's duration)
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DT          (1.0 / 60.0)
#define GRAV        9.8
#define THETA       (15.0 * M_PI / 180.0)
#define MU          0.45
#define RADIUS      0.1
#define V_EPS       0.02
#define REST_SUSTAIN 5
#define MAX_FRAMES  3000
#define JITTER_DEG  3.0
/* contact persistence: a brief airborne gap (chatter from the jittery
 * collider) does not end the contact episode -- only a genuine separation
 * longer than GRACE_FRAMES makes the next touch a "new contact" */
#define GRACE_FRAMES 3
/* contact margin: sub-millimetre chatter from the tilted-normal redirect
 * (v_t * sin(delta) per frame) stays "in contact"; only a real launch
 * separates the body from the surface */
#define SLOP        0.0005

enum { MODE_RAW, MODE_AGE, MODE_GLOBAL };
static const char* mode_name(int m) {
    return m == MODE_RAW ? "raw" : m == MODE_AGE ? "age" : "global";
}

/* per-contact-surface smoothing state (one surface in this scene) */
#define MAX_SURF 8
typedef struct { int used; int age; double nx, ny; } Surf;

static double true_nx, true_ny;   /* surface normal  ( sin15, cos15) */
static double dir_x,  dir_y;      /* surface tangent ( cos15,-sin15) */

/* deterministic per-frame pseudo-random in [0,1): the collider jitter */
static double jitter_hash(int f) {
    double s = sin((double)f * 12.9898) * 43758.5453;
    return s - floor(s);
}

typedef struct {
    int    contact_frames;  /* total frames in contact (contact duration)  */
    int    first_contact;   /* frame index of first contact                */
    int    rest_frames;     /* frames from first contact to rest, -1 never */
    double dangle_deg;      /* normal angle variation (total, degrees)     */
    double slip;            /* slip distance along the surface (m)         */
    double alpha_first;     /* smoothing coefficient at first contact      */
    double alpha_last;      /* smoothing coefficient at last contact frame */
} Metrics;

static Metrics run(int mode, double v0) {
    Metrics m; memset(&m, 0, sizeof m);
    m.first_contact = -1;
    m.rest_frames = -1;
    Surf surfs[MAX_SURF]; memset(surfs, 0, sizeof surfs);

    /* start 0.3 m above the surface at along-slope coordinate 1.0 */
    double px = 1.0 * dir_x + true_nx * (RADIUS + 0.3);
    double py = 1.0 * dir_y + true_ny * (RADIUS + 0.3);
    double vx = v0 * dir_x, vy = v0 * dir_y;

    int episode = 0, gap = 0, rest_run = 0, have_prev = 0;
    double prev_ang = 0.0;

    for (int f = 0; f < MAX_FRAMES; f++) {
        vy -= GRAV * DT;

        /* raw contact normal for this frame (jittered collider) */
        double jr  = (2.0 * jitter_hash(f) - 1.0) * JITTER_DEG * M_PI / 180.0;
        double rnx = true_nx * cos(jr) - true_ny * sin(jr);
        double rny = true_nx * sin(jr) + true_ny * cos(jr);

        double d = px * true_nx + py * true_ny;   /* line through origin */
        int contact = d <= RADIUS + SLOP;

        if (contact) {
            if (m.first_contact < 0) m.first_contact = f;

            Surf* s = &surfs[0];
            if (!episode) {                 /* new contact: reset, no smoothing */
                s->used = 1; s->age = 0; s->nx = rnx; s->ny = rny;
                episode = 1;
            } else {
                s->age++;
            }
            gap = 0;

            double alpha = 1.0;
            if (mode == MODE_AGE)         alpha = 1.0 / (1.0 + s->age);
            else if (mode == MODE_GLOBAL) alpha = 1.0 / (1.0 + f);

            if (mode != MODE_RAW) {
                double nx = s->nx + alpha * (rnx - s->nx);
                double ny = s->ny + alpha * (rny - s->ny);
                double l = hypot(nx, ny);
                if (l > 0) { s->nx = nx / l; s->ny = ny / l; }
            }
            if (m.contact_frames == 0) m.alpha_first = alpha;
            m.alpha_last = alpha;

            double unx = (mode == MODE_RAW) ? rnx : s->nx;
            double uny = (mode == MODE_RAW) ? rny : s->ny;

            double ang = atan2(uny, unx);
            if (have_prev) m.dangle_deg += fabs(ang - prev_ang) * 180.0 / M_PI;
            prev_ang = ang; have_prev = 1;

            /* collision response against the (possibly smoothed) normal */
            double vn = vx * unx + vy * uny;
            if (vn < 0) {
                vx -= unx * vn; vy -= uny * vn;
                /* sliding friction as a force capped by the support load:
                 * mu*N*dt with N = m*g*cos(theta). Jitter impact impulses do
                 * NOT add braking; their real effect is contact intermittency
                 * (the body is kicked off the surface, losing friction). */
                double vt  = vx * dir_x + vy * dir_y;
                double jt  = MU * GRAV * cos(THETA) * DT;
                if (jt > fabs(vt)) jt = fabs(vt);
                double sgn = (vt > 0) - (vt < 0);
                vx -= dir_x * jt * sgn; vy -= dir_y * jt * sgn;
            }
            /* positional projection along the true normal */
            if (d < RADIUS) { px += true_nx * (RADIUS - d); py += true_ny * (RADIUS - d); }

            m.contact_frames++;
            m.slip += fabs(vx * dir_x + vy * dir_y) * DT;

            if (hypot(vx, vy) < V_EPS) {
                if (++rest_run >= REST_SUSTAIN && m.rest_frames < 0)
                    m.rest_frames = f - m.first_contact + 1;
            } else rest_run = 0;
        } else {
            /* brief chatter keeps the episode (and its smoothing state)
             * alive; only a real separation ends it */
            if (++gap > GRACE_FRAMES) { episode = 0; have_prev = 0; }
            rest_run = 0;
        }

        px += vx * DT; py += vy * DT;
        if (m.rest_frames >= 0) break;
        if (px > 100.0 || py < -100.0) break;
    }
    return m;
}

int main(void) {
    true_nx = sin(THETA); true_ny = cos(THETA);
    dir_x   = cos(THETA); dir_y   = -sin(THETA);

    const char*  tier_name[3] = { "short", "medium", "long" };
    const double tier_v0[3]   = { 0.8, 2.0, 4.0 };

    Metrics res[3][3];
    for (int t = 0; t < 3; t++)
        for (int mode = 0; mode < 3; mode++)
            res[t][mode] = run(mode, tier_v0[t]);

    for (int t = 0; t < 3; t++)
        for (int mode = 0; mode < 3; mode++) {
            Metrics* m = &res[t][mode];
            printf("tier=%s mode=%s contact_frames=%d dangle_deg=%.2f slip_m=%.4f rest_frames=%d\n",
                   tier_name[t], mode_name(mode), m->contact_frames,
                   m->dangle_deg, m->slip, m->rest_frames);
        }

    printf("--- nine numbers (mode=age: coefficient from contact duration) ---\n");
    for (int t = 0; t < 3; t++) {
        Metrics* m = &res[t][MODE_AGE];
        printf("%s: dangle_deg=%.2f slip_m=%.4f rest_frames=%d\n",
               tier_name[t], m->dangle_deg, m->slip, m->rest_frames);
    }

    printf("--- baseline (mode=raw) slip for comparison ---\n");
    for (int t = 0; t < 3; t++)
        printf("%s: slip_m=%.4f\n", tier_name[t], res[t][MODE_RAW].slip);

    int reduced = res[0][MODE_AGE].slip < res[0][MODE_RAW].slip;
    printf("short_tier_slip_reduced=%d\n", reduced);
    if (!reduced) {
        /* short-contact slip did not drop: the coefficient was probably
         * derived from global time instead of the contact's own duration */
        printf("diagnostic: coefficient vs contact duration per tier ---\n");
        for (int t = 0; t < 3; t++) {
            Metrics* a = &res[t][MODE_AGE];
            Metrics* g = &res[t][MODE_GLOBAL];
            printf("%s: contact_frames=%d alpha_age[first=%.4f last=%.4f] alpha_global[first=%.4f last=%.4f] slip_global=%.4f\n",
                   tier_name[t], a->contact_frames,
                   a->alpha_first, a->alpha_last,
                   g->alpha_first, g->alpha_last, g->slip);
        }
    }
    return 0;
}
