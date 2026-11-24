
/* bounce: collision response at an inside corner.
 *
 * Two constraints keep the rebound physically self-consistent:
 *   1. restitution is derived from the incoming velocity and the contact
 *      normal (normal component), not from a fixed coefficient;
 *   2. the rebound normal speed is clamped to the incoming normal speed,
 *      and when the clamp bites the velocity is rebuilt along the reflected
 *      direction instead of just dropping the normal component.
 * A contact shorter than MIN_CONTACT frames only gets positional projection
 * (velocity untouched), never a full rebound. */
#include <stdio.h>
#include <math.h>

#define PI_F 3.14159265358979323846f
#define DT 0.05f
#define MIN_CONTACT 2
#define MAX_FRAMES 200
#define STUCK_LIMIT 30
#define GATE_VN 1e-4f

/* restitution model: grows with impact speed, not bounded by construction */
#define REST_E0 0.40f
#define REST_K  0.045f

/* inside corner: two perpendicular walls meeting at the origin */
#define CONTACT_EPS 1e-4f
#define START_GAP 0.02f
#define DIR_N (-0.948683298050514f)  /* inward dir, normal component: -3/sqrt10 */
#define DIR_T (-0.316227766016838f)  /* tangential component:       -1/sqrt10 */

typedef struct { float x, y; } V2;
typedef enum { POLICY_FIXED = 0, POLICY_NAIVE = 1 } Policy;

static float dot2(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static float len2(V2 a) { return sqrtf(dot2(a, a)); }
static V2 scale2(V2 a, float s) { return (V2){ a.x * s, a.y * s }; }
static V2 add2(V2 a, V2 b) { return (V2){ a.x + b.x, a.y + b.y }; }
static V2 sub2(V2 a, V2 b) { return (V2){ a.x - b.x, a.y - b.y }; }

/* corner contact normal: bisector when both walls are penetrated,
 * otherwise the normal of the single wall that is penetrated */
static int corner_normal(V2 p, V2* n) {
    int hx = p.x < CONTACT_EPS;
    int hy = p.y < CONTACT_EPS;
    if (!hx && !hy) return 0;
    if (hx && hy) *n = (V2){ 0.707106781186548f, 0.707106781186548f };
    else if (hx)  *n = (V2){ 1.0f, 0.0f };
    else          *n = (V2){ 0.0f, 1.0f };
    return 1;
}

/* positional projection: push out of whichever wall is actually penetrated.
 * independent of the (locked) response normal, so sliding along one wall
 * after a corner rebound does not get projected along the corner bisector */
static void project_out(V2* p) {
    if (p->x < CONTACT_EPS) p->x = CONTACT_EPS;
    if (p->y < CONTACT_EPS) p->y = CONTACT_EPS;
}

/* velocity-dependent restitution from incoming speed along the normal */
static float restitution_raw(float vn_in) {
    return REST_E0 + REST_K * fabsf(vn_in);
}

/* full rebound before the speed cap; shared by both policies so the
 * pre-clamp quantities in the diagnostic table are identical */
static void respond_raw(V2 v, V2 n, V2* v_pre, float* vn_in, float* vn_reb) {
    float vn = dot2(v, n);
    V2 vt = sub2(v, scale2(n, vn));
    float reb = -restitution_raw(vn) * vn;
    *vn_in = vn;
    *vn_reb = reb;
    *v_pre = add2(vt, scale2(n, reb));
}

/* FIXED: clamp the rebound normal speed, then rebuild along the reflection */
static V2 respond_fixed(V2 v, V2 n) {
    V2 v_pre; float vn_in, vn_reb;
    respond_raw(v, n, &v_pre, &vn_in, &vn_reb);
    float cap = fabsf(vn_in);
    if (vn_reb > cap) vn_reb = cap;
    V2 vt = sub2(v, scale2(n, vn_in));
    return add2(vt, scale2(n, vn_reb));
}

/* NAIVE (buggy baseline): when the cap bites, the normal component is
 * discarded rather than reflected, so the body keeps re-penetrating */
static V2 respond_naive(V2 v, V2 n) {
    V2 v_pre; float vn_in, vn_reb;
    respond_raw(v, n, &v_pre, &vn_in, &vn_reb);
    if (vn_reb <= fabsf(vn_in)) return v_pre;
    return sub2(v, scale2(n, vn_in));
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("check %s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

typedef struct {
    float in_speed, ratio, angle_deg;
    int stuck_frames;
    V2 n, v_in, v_pre, v_post;
} Result;

static Result run_impact(float speed, Policy policy) {
    Result r = { 0 };
    V2 p = { START_GAP, START_GAP };
    V2 v = { DIR_N * speed, DIR_T * speed };
    V2 n = { 0 };
    int contact = 0, measured = 0, first = 1, stuck = 0;

    r.in_speed = speed;
    for (int frame = 0; frame < MAX_FRAMES && stuck < STUCK_LIMIT; frame++) {
        V2 cn;
        if (corner_normal(p, &cn)) {
            if (!contact) n = cn; /* normal fixed for the whole episode */
            contact++;
            float vn = dot2(v, n);

            if (vn < -GATE_VN) {
                if (first) { /* pristine incoming state, pre-clamp response */
                    float vn_in2, vn_reb2;
                    respond_raw(v, n, &r.v_pre, &vn_in2, &vn_reb2);
                    r.v_in = v;
                    first = 0;
                }
                if (contact >= MIN_CONTACT) {
                    v = (policy == POLICY_FIXED) ? respond_fixed(v, n)
                                                : respond_naive(v, n);
                    if (!measured) {
                        r.v_post = v;
                        r.ratio = len2(v) / speed;
                        float c = dot2(v, n) / len2(v);
                        if (c > 1.0f) c = 1.0f;
                        r.angle_deg = acosf(c) * 180.0f / PI_F;
                        measured = 1;
                    }
                }
                /* gate still closed: project out of the wall, keep velocity */
            } else if (contact >= MIN_CONTACT && fabsf(vn) < GATE_VN) {
                /* pressed against the wall with no normal motion at all */
                stuck++;
            }

            project_out(&p);
        } else {
            contact = 0;
            if (measured || !first) break; /* episode over */
        }
        p = add2(p, scale2(v, DT));
    }

    if (!measured) { /* never escaped: pre-clamp response is the diagnostic */
        r.v_post = r.v_pre;
        r.ratio = len2(r.v_pre) / speed;
        float c = dot2(r.v_pre, n) / len2(r.v_pre);
        if (c > 1.0f) c = 1.0f;
        r.angle_deg = acosf(c) * 180.0f / PI_F;
    }
    r.stuck_frames = stuck;
    r.n = n;
    return r;
}

static void print_table(const char* label, const Result* rs, int n) {
    printf("policy=%s\n", label);
    for (int i = 0; i < n; i++) {
        printf("  v_in=%6.2f  ratio=%.4f  stuck=%3d  angle_deg=%6.2f\n",
               rs[i].in_speed, rs[i].ratio, rs[i].stuck_frames,
               rs[i].angle_deg);
    }
}

int main(void) {
    /* ---- unit checks ---- */
    V2 n = { 0.707106781186548f, 0.707106781186548f };

    check("restitution_speed_dependent",
          fabsf(restitution_raw(2.0f) - restitution_raw(10.0f)) > 1e-6f);

    V2 vi = { -3.0f, -1.0f };
    V2 vf = respond_fixed(vi, n);
    check("fixed_reflects_along_normal", dot2(vf, n) > 0.0f);
    check("fixed_rebound_not_faster", len2(vf) <= len2(vi) + 1e-6f);

    /* high-speed impact whose raw restitution exceeds 1.0 -> cap bites */
    V2 hi = { -18.0f, -6.0f };
    V2 hp; float vn_in, vn_reb;
    respond_raw(hi, n, &hp, &vn_in, &vn_reb);
    check("cap_engages_at_high_speed", vn_reb > fabsf(vn_in));
    V2 hf = respond_fixed(hi, n);
    check("fixed_clamp_keeps_direction", dot2(hf, n) > 0.0f);
    /* pure normal incidence: naive response leaves no normal component */
    V2 hn_in = scale2(n, -20.0f);
    V2 hnaive = respond_naive(hn_in, n);
    check("naive_cap_kills_normal", fabsf(dot2(hnaive, n)) < 1e-6f);

    /* contact gate: a contact ending before MIN_CONTACT never alters velocity */
    {
        V2 q = { CONTACT_EPS * 0.5f, START_GAP };
        V2 w = { 0.5f, 1.0f }, w0 = w;
        V2 wn; int touched = 0;
        for (int k = 0; k < MIN_CONTACT - 1; k++) {
            if (corner_normal(q, &wn)) {
                touched++;
                project_out(&q);
            }
            q = add2(q, scale2(w, DT));
        }
        check("gate_holds_rebound",
              touched == 1 && w.x == w0.x && w.y == w0.y);
    }

    /* ---- three-speed experiment, nine numbers per policy ---- */
    float speeds[3] = { 2.0f, 6.0f, 18.0f };
    Result fixed[3], naive[3];
    for (int i = 0; i < 3; i++) {
        fixed[i] = run_impact(speeds[i], POLICY_FIXED);
        naive[i] = run_impact(speeds[i], POLICY_NAIVE);
    }

    print_table("fixed", fixed, 3);
    print_table("naive", naive, 3);

    /* ---- diagnostic: high-speed tier stuck most => cap without direction fix */
    int hi_stuck_worst =
        naive[2].stuck_frames > naive[0].stuck_frames &&
        naive[2].stuck_frames >= naive[1].stuck_frames;
    printf("diagnostic high_speed_stuck_worst=%d\n", hi_stuck_worst ? 1 : 0);
    if (hi_stuck_worst) {
        printf("clamp_table (naive policy, at first contact frame):\n");
        printf("  tier   v_in    vn_in   e_raw   vn_pre     cap   vn_post   "
               "vx_pre   vy_pre   vx_post  vy_post\n");
        for (int i = 0; i < 3; i++) {
            Result* r = &naive[i];
            float vn_in2 = dot2(r->v_in, r->n);
            float e = restitution_raw(vn_in2);
            float vn_pre = -e * vn_in2;
            float cap = fabsf(vn_in2);
            float vn_post = dot2(r->v_post, r->n);
            printf("  %s %6.2f %8.4f %6.3f %9.4f %8.4f %9.4f %8.4f %8.4f "
                   "%8.4f %8.4f\n",
                   i == 0 ? "low " : i == 1 ? "mid " : "high",
                   r->in_speed, vn_in2, e, vn_pre, cap, vn_post,
                   r->v_pre.x, r->v_pre.y, r->v_post.x, r->v_post.y);
        }
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
