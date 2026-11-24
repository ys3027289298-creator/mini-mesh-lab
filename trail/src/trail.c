
/* trail: swing-trail anchor point from support functions.
 *
 * Support function of shape S along unit direction d:
 *     h_S(d) = max_{p in S} <p, d>,   s_S(d) = argmax
 * Vertex mesh: max over vertices (NOT bbox center + radius). Rounded
 * (arc) edge: max over the arc, which is generally not a vertex.
 *
 * Trail start p0 = s_S(axis). One swing rotates the weapon THETA rad
 * about the pivot; the trail is the head-point polyline capped at
 * BUDGET world units:
 *     r = |p0 - O|,  L(p0) = min(r*THETA, BUDGET)
 * Moving p0 outward raises r, so a capped budget covers angle B/r:
 * anchor accuracy and coverage pull against each other. Settle both
 * together: p0 from the (face-aware) support function, budget sized
 * as B >= r*THETA for full coverage.
 */
#include <stdio.h>
#include <math.h>

#define MAXV 16
#define THETA  2.0943951f   /* swing angle, 120 deg */
#define BUDGET 2.5f         /* trail length budget, world units */

typedef struct { float x, y, z; } V3;

static V3 add(V3 p, V3 q) { V3 r = { p.x+q.x, p.y+q.y, p.z+q.z }; return r; }
static V3 sub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 scale(V3 p, float s) { V3 r = { p.x*s, p.y*s, p.z*s }; return r; }
static V3 cross(V3 p, V3 q) { V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r; }
static float dot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float len(V3 p) { return sqrtf(dot(p, p)); }
static V3 norm(V3 p) { float l = len(p); V3 r = {0,0,0}; if (l > 0) { r.x=p.x/l; r.y=p.y/l; r.z=p.z/l; } return r; }

typedef struct {
    V3 c; float r;   /* arc center, radius */
    V3 n;            /* plane normal (unit) */
    V3 m;            /* mid-arc direction in plane (unit) */
    float half;      /* half-angle of arc span, rad */
    int valid;
} Arc;

typedef struct {
    const char* name;
    V3 v[MAXV]; int nv;
    V3 axis;         /* authored main axis (unit) */
    V3 pivot;        /* swing pivot */
    Arc edge;        /* rounded cutting edge; valid=0 for polygonal tips */
} Weapon;

/* vertex-only support: farthest mesh vertex along d */
static V3 support_vertices(const Weapon* w, V3 d, float* h) {
    int bi = 0; float bh = -1e30f;
    for (int i = 0; i < w->nv; i++) {
        float s = dot(w->v[i], d);
        if (s > bh) { bh = s; bi = i; }
    }
    *h = bh;
    return w->v[bi];
}

/* support of a circular arc edge: c + r*u with u the in-plane direction
 * of d, clamped to the arc span -- generally lies between vertices */
static V3 support_arc(const Arc* a, V3 d, float* h) {
    V3 dp = sub(d, scale(a->n, dot(d, a->n)));
    float l = len(dp);
    V3 u = (l < 1e-9f) ? a->m : scale(dp, 1.0f / l);
    V3 t = cross(a->n, a->m);            /* in-plane, 90 deg from m */
    float ca = dot(u, a->m);
    if (ca > 1.0f) ca = 1.0f;
    if (ca < -1.0f) ca = -1.0f;
    if (acosf(ca) > a->half) {           /* outside span: nearer endpoint */
        float side = dot(u, t) >= 0.0f ? 1.0f : -1.0f;
        float phi = side * a->half;
        u = add(scale(a->m, cosf(phi)), scale(t, sinf(phi)));
    }
    V3 p = add(a->c, scale(u, a->r));
    *h = dot(p, d);
    return p;
}

/* face-aware support: best of vertices and rounded edge */
static V3 support_true(const Weapon* w, V3 d, float* h) {
    float hv; V3 pv = support_vertices(w, d, &hv);
    if (!w->edge.valid) { *h = hv; return pv; }
    float ha; V3 pa = support_arc(&w->edge, d, &ha);
    if (ha > hv) { *h = ha; return pa; }
    *h = hv; return pv;
}

static V3 rot_z_about(V3 p, V3 pivot, float a) {
    float c = cosf(a), s = sinf(a);
    V3 q = sub(p, pivot);
    V3 r = { c*q.x - s*q.y + pivot.x, s*q.x + c*q.y + pivot.y, q.z + pivot.z };
    return r;
}

/* one swing: head starts at p0, rotates THETA about pivot; visible
 * length = accumulated head-polyline length, capped at budget */
static float swing_visible_len(V3 p0, V3 pivot, float theta, float budget) {
    const int steps = 256;
    V3 prev = p0;
    float acc = 0.0f;
    for (int i = 1; i <= steps; i++) {
        V3 p = rot_z_about(p0, pivot, theta * (float)i / (float)steps);
        float seg = len(sub(p, prev));
        if (acc + seg >= budget) return budget;
        acc += seg;
        prev = p;
    }
    return acc;
}

/* bbox center + radius along axis: the estimate we are replacing */
static float center_radius_est(const Weapon* w, V3 d) {
    V3 mn = w->v[0], mx = w->v[0];
    for (int i = 1; i < w->nv; i++) {
        V3 p = w->v[i];
        if (mn.x > p.x) mn.x = p.x;
        if (mn.y > p.y) mn.y = p.y;
        if (mn.z > p.z) mn.z = p.z;
        if (mx.x < p.x) mx.x = p.x;
        if (mx.y < p.y) mx.y = p.y;
        if (mx.z < p.z) mx.z = p.z;
    }
    V3 c = scale(add(mn, mx), 0.5f);
    return dot(c, d) + 0.5f * len(sub(mx, mn));
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Weapon ws[3];

    /* straight sword: polygonal, tip vertex on axis */
    ws[0].name = "straight_sword";
    ws[0].v[0] = (V3){-0.03f, 0.20f, 0}; ws[0].v[1] = (V3){0.03f, 0.20f, 0};
    ws[0].v[2] = (V3){-0.03f, 0.85f, 0}; ws[0].v[3] = (V3){0.03f, 0.85f, 0};
    ws[0].v[4] = (V3){0.00f, 1.00f, 0};
    ws[0].nv = 5;
    ws[0].axis = (V3){0, 1, 0};
    ws[0].pivot = (V3){0, 0, 0};
    ws[0].edge.valid = 0;

    /* curved blade: polygonal curve, tip vertex off the world axes */
    ws[1].name = "curved_blade";
    ws[1].v[0] = (V3){-0.025f, 0.20f, 0}; ws[1].v[1] = (V3){0.025f, 0.20f, 0};
    ws[1].v[2] = (V3){0.00f, 0.50f, 0};  ws[1].v[3] = (V3){0.08f, 0.72f, 0};
    ws[1].v[4] = (V3){0.15f, 0.85f, 0};  ws[1].v[5] = (V3){0.25f, 0.93f, 0};
    ws[1].nv = 6;
    ws[1].axis = norm((V3){0.25f, 0.93f, 0});
    ws[1].pivot = (V3){0, 0, 0};
    ws[1].edge.valid = 0;

    /* wide blade: rounded edge is a true circular arc face, r=0.35,
     * half-span 60 deg; only its two endpoints are mesh vertices */
    ws[2].name = "wide_blade";
    ws[2].v[0] = (V3){-0.15f, 0.20f, 0};      ws[2].v[1] = (V3){0.15f, 0.20f, 0};
    ws[2].v[2] = (V3){-0.28f, 0.60f, 0};      ws[2].v[3] = (V3){0.28f, 0.60f, 0};
    ws[2].v[4] = (V3){-0.3031089f, 0.925f, 0}; ws[2].v[5] = (V3){0.3031089f, 0.925f, 0};
    ws[2].nv = 6;
    ws[2].axis = (V3){0, 1, 0};
    ws[2].pivot = (V3){0, 0, 0};
    ws[2].edge.c = (V3){0, 0.75f, 0};
    ws[2].edge.r = 0.35f;
    ws[2].edge.n = (V3){0, 0, 1};
    ws[2].edge.m = (V3){0, 1, 0};
    ws[2].edge.half = 1.0471976f;   /* 60 deg */
    ws[2].edge.valid = 1;

    float dev[3], tlen[3], fit[3];
    V3 pv[3], pt[3];
    float hv[3], ht[3];

    for (int k = 0; k < 3; k++) {
        Weapon* w = &ws[k];
        pv[k] = support_vertices(w, w->axis, &hv[k]);  /* as implemented: vertex-only */
        pt[k] = support_true(w, w->axis, &ht[k]);      /* ground truth: face-aware */
        float base = 1e30f;
        for (int i = 0; i < w->nv; i++) {
            float s = dot(w->v[i], w->axis);
            if (base > s) base = s;
        }
        dev[k] = ht[k] - hv[k];                        /* axial anchor error */
        tlen[k] = swing_visible_len(pv[k], w->pivot, THETA, BUDGET);
        float blade = ht[k] - base;
        float off = len(sub(pv[k], pt[k]));            /* anchor offset from true tip */
        fit[k] = blade > 0.0f ? 1.0f - off / blade : 1.0f;
        if (fit[k] < 0.0f) fit[k] = 0.0f;

        printf("[%s]\n", w->name);
        printf("start_dev=%.4f\n", dev[k]);
        printf("trail_len=%.4f\n", tlen[k]);
        printf("fit=%.4f\n", fit[k]);
        printf("center_radius_est=%.4f\n", center_radius_est(w, w->axis));
    }

    /* if the wide blade fits worst, vertex-only support missed the arc
     * edge face: put vertex farthest and face farthest side by side */
    int worst = 0;
    for (int k = 1; k < 3; k++) if (fit[k] < fit[worst]) worst = k;
    if (ws[worst].edge.valid) {
        printf("[%s] worst_fit: vertex-only support missed the arc edge\n", ws[worst].name);
        printf("vertex_farthest=(%.4f,%.4f,%.4f) axial=%.4f\n",
               pv[worst].x, pv[worst].y, pv[worst].z, hv[worst]);
        printf("face_farthest=(%.4f,%.4f,%.4f) axial=%.4f\n",
               pt[worst].x, pt[worst].y, pt[worst].z, ht[worst]);
    }

    check("straight_dev_zero", fabsf(dev[0]) < 1e-6f);
    check("curved_dev_zero", fabsf(dev[1]) < 1e-6f);
    check("wide_dev_arc", fabsf(dev[2] - 0.175f) < 1e-3f);
    check("len_formula_straight", fabsf(tlen[0] - 1.0f * THETA) < 1e-3f);
    check("len_formula_wide", fabsf(tlen[2] - len(pv[2]) * THETA) < 1e-3f);
    check("wide_worst_fit", fit[2] < fit[0] && fit[2] < fit[1]);
    check("budget_not_binding", tlen[0] < BUDGET && tlen[1] < BUDGET && tlen[2] < BUDGET);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
