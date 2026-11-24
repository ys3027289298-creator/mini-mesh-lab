/* landing: ragdoll ground landing.
 *
 * The landing pose is solved in two fixed stages:
 *   1. align_normal - rotate the body about the root origin so its local +Z
 *      (the chest normal of a supine body) matches the terrain normal.
 *   2. clamp_root   - translate the root along world +Z so the lowest point
 *      of the root bone has zero clearance to the terrain.
 *
 * The stage order is enforced by run_landing() and recorded in Trace.  The
 * penetration metrics are measurements of the result; they are not used as a
 * proxy for the order.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846
#define MAXP 32
#define ROOT_RADIUS 0.12   /* pelvis half-thickness: root lowest point is at -r */
#define PEN_EPS 1e-4       /* m; deeper than this counts as penetrating */
#define CONTACT_EPS 1e-3   /* m; within this band counts as touching */

typedef struct { double x, y, z; } V3;

/* z = tan(deg) * x, a plane rising toward +X. */
typedef struct { double deg, sin_t, cos_t; } Terrain;

typedef struct {
    V3 local[MAXP];  /* relative to the root origin; [0] is the root lowest point */
    int n;
    V3 root;         /* world position of the root origin */
    double angle;    /* rotation about world Y; +angle maps +Z to the uphill normal */
} Pose;

typedef enum { STAGE_ALIGN_NORMAL = 1, STAGE_CLAMP_ROOT = 2 } Stage;

typedef struct { Stage ev[4]; int n; } Trace;

typedef struct {
    int penetrating;      /* vertices with clearance < -PEN_EPS */
    double max_depth;     /* max(0, -clearance), metres */
    int contacts;         /* vertices with |clearance| <= CONTACT_EPS */
    double min_clearance; /* signed, metres */
} Metrics;

static V3 sub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 cross(V3 p, V3 q) { V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r; }
static double dot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static double len(V3 p) { return sqrt(dot(p,p)); }

static Terrain terrain_make(double deg) {
    Terrain t;
    double r = deg * PI / 180.0;
    t.deg = deg;
    t.sin_t = sin(r);
    t.cos_t = cos(r);
    return t;
}

static double terrain_height(const Terrain* t, double x, double y) {
    (void)y;
    return (t->sin_t / t->cos_t) * x;
}

static V3 terrain_normal(const Terrain* t) {
    V3 n = { -t->sin_t, 0.0, t->cos_t };
    return n;
}

/* Signed distance from the terrain plane along the unit normal. */
static double clearance(const Terrain* t, V3 p) {
    return (p.z - terrain_height(t, p.x, p.y)) * t->cos_t;
}

static int pose_add(Pose* p, double x, double y, double z) {
    if (p->n >= MAXP) return -1;
    p->local[p->n].x = x;
    p->local[p->n].y = y;
    p->local[p->n].z = z;
    return p->n++;
}

/* Supine body, +X head, +Y left, +Z chest normal.  The back is a flat
 * support plane at local z = -ROOT_RADIUS, so the root lowest point is also
 * the lowest point of the body in the landing pose. */
static void body_init(Pose* p) {
    memset(p, 0, sizeof(*p));
    pose_add(p, 0.00,  0.00, -ROOT_RADIUS);  /* 0: root bone lowest point */
    pose_add(p, 0.62,  0.00, -ROOT_RADIUS);  /* 1: back of head */
    pose_add(p, 0.38,  0.19, -ROOT_RADIUS);  /* 2: shoulder L */
    pose_add(p, 0.38, -0.19, -ROOT_RADIUS);  /* 3: shoulder R */
    pose_add(p, 0.10,  0.00, -ROOT_RADIUS);  /* 4: mid back */
    pose_add(p, 0.00,  0.13, -ROOT_RADIUS);  /* 5: pelvis L */
    pose_add(p, 0.00, -0.13, -ROOT_RADIUS);  /* 6: pelvis R */
    pose_add(p, -0.42,  0.09, -ROOT_RADIUS); /* 7: knee L */
    pose_add(p, -0.42, -0.09, -ROOT_RADIUS); /* 8: knee R */
    pose_add(p, -0.80,  0.07, -ROOT_RADIUS); /* 9: heel L */
    pose_add(p, -0.80, -0.07, -ROOT_RADIUS); /* 10: heel R */

    pose_add(p, 0.62,  0.00,  0.06);         /* 11: crown */
    pose_add(p, 0.30,  0.00,  0.08);         /* 12: chest */
    pose_add(p, 0.05,  0.00,  0.05);         /* 13: abdomen */
    pose_add(p, 0.20,  0.38, -0.02);         /* 14: hand L */
    pose_add(p, 0.20, -0.38, -0.02);         /* 15: hand R */
    pose_add(p, -0.84,  0.07, -0.02);        /* 16: toe L */
    pose_add(p, -0.84, -0.07, -0.02);        /* 17: toe R */
}

static V3 pose_point(const Pose* p, int i) {
    const V3 l = p->local[i];
    double ca = cos(p->angle);
    double sa = sin(p->angle);
    V3 w;
    w.x = l.x * ca - l.z * sa + p->root.x;
    w.y = l.y + p->root.y;
    w.z = l.x * sa + l.z * ca + p->root.z;
    return w;
}

static void trace_reset(Trace* tr) { tr->n = 0; }

static void trace_push(Trace* tr, Stage s) {
    if (tr->n < (int)(sizeof(tr->ev) / sizeof(tr->ev[0]))) tr->ev[tr->n++] = s;
}

static int trace_is(const Trace* tr, Stage first, Stage second) {
    return tr->n == 2 && tr->ev[0] == first && tr->ev[1] == second;
}

static const char* trace_name(const Trace* tr) {
    if (tr->n != 2) return "incomplete";
    if (tr->ev[0] == STAGE_ALIGN_NORMAL && tr->ev[1] == STAGE_CLAMP_ROOT)
        return "align_normal->clamp_root";
    if (tr->ev[0] == STAGE_CLAMP_ROOT && tr->ev[1] == STAGE_ALIGN_NORMAL)
        return "clamp_root->align_normal";
    return "other";
}

/* Stage 1: rotate about the root origin so local +Z matches the terrain
 * normal at the landing point.  No translation in this stage. */
static void stage_align_normal(Pose* p, const Terrain* t, Trace* tr) {
    p->angle = t->deg * PI / 180.0;
    trace_push(tr, STAGE_ALIGN_NORMAL);
}

/* Stage 2: move the root along world +Z until the root bone's lowest point
 * has zero clearance.  A world +Z move of dz changes clearance by dz*cos_t. */
static void stage_clamp_root(Pose* p, const Terrain* t, Trace* tr) {
    double c = clearance(t, pose_point(p, 0));
    p->root.z += -c / t->cos_t;
    trace_push(tr, STAGE_CLAMP_ROOT);
}

static Metrics measure(const Pose* p, const Terrain* t) {
    Metrics m;
    m.penetrating = 0;
    m.max_depth = 0.0;
    m.contacts = 0;
    m.min_clearance = 1e30;

    for (int i = 0; i < p->n; i++) {
        double c = clearance(t, pose_point(p, i));
        if (c < m.min_clearance) m.min_clearance = c;
        if (c < -PEN_EPS) m.penetrating++;
        if (-c > m.max_depth) m.max_depth = -c;
        if (fabs(c) <= CONTACT_EPS) m.contacts++;
    }
    if (m.max_depth < 0.0) m.max_depth = 0.0;
    return m;
}

/* reverse_control exists only for the order-probe self-check; the shipped
 * landing path always uses align then clamp. */
static void run_landing(double deg, int reverse_control, Pose* out, Metrics* m, Trace* tr) {
    Terrain t = terrain_make(deg);
    Pose p;
    body_init(&p);
    trace_reset(tr);

    if (reverse_control) {
        stage_clamp_root(&p, &t, tr);
        stage_align_normal(&p, &t, tr);
    } else {
        stage_align_normal(&p, &t, tr);
        stage_clamp_root(&p, &t, tr);
    }

    *out = p;
    *m = measure(&p, &t);
}

static const double SLOPES[3] = { 0.0, 15.0, 35.0 };

static void print_table(const char* title, int reverse,
                        Pose ps[3], Metrics mm[3], Trace tr[3]) {
    printf("%s\n", title);
    printf("%-10s %-21s %-19s %-15s %s\n",
           "slope_deg", "penetrating_vertices", "max_penetration_m",
           "contact_points", "actual_order");
    for (int i = 0; i < 3; i++) {
        run_landing(SLOPES[i], reverse, &ps[i], &mm[i], &tr[i]);
        printf("%-10.0f %-21d %-19.6f %-15d %s\n",
               SLOPES[i], mm[i].penetrating, mm[i].max_depth,
               mm[i].contacts, trace_name(&tr[i]));
    }
}

static int failures = 0;

static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Pose ps[3], pr[3];
    Metrics mm[3], mr[3];
    Trace tr[3], trr[3];

    print_table("main (required order)", 0, ps, mm, tr);
    printf("\n");
    print_table("order probe (reverse-order control; NOT the shipped order)",
                1, pr, mr, trr);
    printf("\n");

    int steep_main = mm[2].max_depth > mm[0].max_depth + 1e-12 &&
                     mm[2].max_depth > mm[1].max_depth + 1e-12;
    int steep_rev  = mr[2].max_depth > mr[0].max_depth + 1e-12 &&
                     mr[2].max_depth > mr[1].max_depth + 1e-12;
    printf("main_steep_is_deepest=%d\n", steep_main);
    printf("reverse_steep_is_deepest=%d\n", steep_rev);

    int order_ok = 1, no_pen = 1, contacts_ok = 1, contacts_same = 1;
    int root_on_ground = 1;
    for (int i = 0; i < 3; i++) {
        order_ok &= trace_is(&tr[i], STAGE_ALIGN_NORMAL, STAGE_CLAMP_ROOT);
        no_pen &= (mm[i].penetrating == 0 && mm[i].max_depth <= 1e-9);
        contacts_ok &= (mm[i].contacts > 0);
        contacts_same &= (mm[i].contacts == mm[0].contacts);
        Terrain t = terrain_make(SLOPES[i]);
        root_on_ground &= fabs(clearance(&t, pose_point(&ps[i], 0))) < 1e-9;
    }
    check("main_order_align_then_clamp", order_ok);
    check("main_no_penetration", no_pen);
    check("main_has_contacts", contacts_ok);
    check("main_contacts_slope_invariant", contacts_same);
    check("main_root_low_point_on_ground", root_on_ground);
    check("main_steep_not_deepest", steep_main == 0);

    int probe_order_ok = 1;
    for (int i = 0; i < 3; i++)
        probe_order_ok &= trace_is(&trr[i], STAGE_CLAMP_ROOT, STAGE_ALIGN_NORMAL);
    check("probe_order_clamp_then_align", probe_order_ok);
    check("probe_steep_is_deepest", steep_rev == 1);
    check("probe_depth_monotonic_with_slope",
          mr[2].max_depth > mr[1].max_depth && mr[1].max_depth > mr[0].max_depth);
    check("probe_steep_penetrates", mr[2].penetrating > 0);

    double expect = ROOT_RADIUS * (1.0 - cos(SLOPES[2] * PI / 180.0));
    check("probe_depth_matches_root_offset", fabs(mr[2].max_depth - expect) < 1e-6);

    Terrain t35 = terrain_make(SLOPES[2]);
    V3 n = terrain_normal(&t35);
    check("terrain_normal_unit", fabs(len(n) - 1.0) < 1e-12);
    check("flat_terrain_zero_height", fabs(terrain_height(&t35, 0.0, 0.0)) < 1e-12);
    check("slope_height_matches_tan", fabs(terrain_height(&t35, 1.0, 0.0) -
                                           tan(SLOPES[2] * PI / 180.0)) < 1e-12);

    /* Keep the frame/terrain helpers exercised (silences unused warnings for
     * sub/cross/dot in builds where only the self-check changes). */
    V3 a = { 1.0, 0.0, 0.0 }, b = { 0.0, 1.0, 0.0 };
    (void)sub(a, b);
    (void)cross(a, b);
    (void)dot(a, b);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
