
/*
 * proxy: distance-clustered collision proxies for hundreds of destructibles.
 *
 * Model (deliberately small and fully deterministic):
 *   - 2D top-down collision world; every destructible is a yawed box (OBB).
 *   - objects closer than CLUSTER_R are union-find merged into one proxy.
 *   - each proxy shape is either an AABB or a monotone-chain convex hull.
 *   - hull is chosen only when BOTH gates pass:
 *       concavity low enough : C <= CONCAVE_THR   (concavity_gate)
 *       arrangement tight    : T <= TIGHT_THR     (tightness_gate)
 *     otherwise the proxy falls back to an AABB.
 *   - per-frame workload: BENCH_RAYS collision probes (projectile rays).
 *   - ground truth = rays tested against every individual OBB.
 *
 * AABB and hull are built from all member corners, so with an unlimited
 * narrow phase a proxy can never geometrically miss an OBB hit
 * (check_conservative_proxies). Non-zero "collision misses" therefore model
 * a real engine mechanism: the per-frame narrow-phase budget NARROW_BUDGET
 * is exhausted by loose proxies, so unserved probes are dropped; a dropped
 * probe the oracle proves was a real hit counts as a miss.
 *
 * The benchmark is run with the correct concavity sign and with an inverted
 * control (C >= threshold). The inverted control must reproduce the
 * diagnostic signature "high density has the most misses".
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define MAX_OBJ      640
#define MAX_PROXY    96
#define MAX_MEMBERS  128
#define MAX_PTS      512
#define BENCH_RAYS   512
#define BENCH_FRAMES 120
#define GEO_FRAMES   6
#define TIME_REPS    5
#define TIME_FRAMES  120
#define NARROW_BUDGET 8192L
#define SAMPLE_GRID  28

#define CLUSTER_R    4.20f
#define CONCAVE_THR  0.60f
#define TIGHT_THR    1.80f

typedef struct { float x, y; } V2;

/* ---------------- deterministic RNG (xorshift32) ---------------- */
typedef struct { uint32_t s; } Rng;

static float rng_u(Rng* r) {
    uint32_t x = r->s ? r->s : 1u;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    r->s = x;
    return (x >> 8) * (1.0f / 16777216.0f); /* [0,1), 24-bit mantissa */
}
static float rng_range(Rng* r, float a, float b) { return a + (b - a) * rng_u(r); }

/* ---------------- geometry ---------------- */
static V2 vsub(V2 a, V2 b) { return (V2){ a.x - b.x, a.y - b.y }; }

static float cross2(V2 o, V2 a, V2 b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

typedef struct {
    V2 c;
    float ca, sa;
    float hx, hy;
    V2 corner[4];
} Obj;

static void obj_init(Obj* o, float x, float y, float ang, float hx, float hy) {
    o->c = (V2){ x, y };
    o->ca = cosf(ang);
    o->sa = sinf(ang);
    o->hx = hx;
    o->hy = hy;
    static const float sx[4] = { -1, 1, 1, -1 };
    static const float sy[4] = { -1, -1, 1, 1 };
    for (int i = 0; i < 4; i++) {
        float lx = sx[i] * hx, ly = sy[i] * hy;
        o->corner[i] = (V2){ x + lx * o->ca - ly * o->sa,
                             y + lx * o->sa + ly * o->ca };
    }
}

static int ray_aabb(V2 o, V2 d, V2 mn, V2 mx) {
    float tmin = -1e30f, tmax = 1e30f;
    float ox[2] = { o.x, o.y }, dx[2] = { d.x, d.y };
    float lo[2] = { mn.x, mn.y }, hi[2] = { mx.x, mx.y };
    for (int i = 0; i < 2; i++) {
        if (fabsf(dx[i]) < 1e-12f) {
            if (ox[i] < lo[i] || ox[i] > hi[i]) return 0;
        } else {
            float t1 = (lo[i] - ox[i]) / dx[i];
            float t2 = (hi[i] - ox[i]) / dx[i];
            if (t1 > t2) { float t = t1; t1 = t2; t2 = t; }
            if (t1 > tmin) tmin = t1;
            if (t2 < tmax) tmax = t2;
            if (tmin > tmax) return 0;
        }
    }
    return tmax >= 0.0f && tmin <= 1e9f;
}

static int ray_obb(V2 ro, V2 rd, const Obj* o) {
    V2 ol = { (ro.x - o->c.x) * o->ca + (ro.y - o->c.y) * o->sa,
             -(ro.x - o->c.x) * o->sa + (ro.y - o->c.y) * o->ca };
    V2 dl = { rd.x * o->ca + rd.y * o->sa,
             -rd.x * o->sa + rd.y * o->ca };
    return ray_aabb(ol, dl, (V2){ -o->hx, -o->hy }, (V2){ o->hx, o->hy });
}

static int point_in_obb(V2 p, const Obj* o) {
    float dx = p.x - o->c.x, dy = p.y - o->c.y;
    float lx = dx * o->ca + dy * o->sa;
    float ly = -dx * o->sa + dy * o->ca;
    return fabsf(lx) <= o->hx && fabsf(ly) <= o->hy;
}

static int cmp_v2(const void* pa, const void* pb) {
    const V2* a = pa, * b = pb;
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    return 0;
}

/* monotone chain; returns hull vertex count (ccw) */
static int hull2d(const V2* pts, int n, V2* out) {
    if (n > MAX_PTS) n = MAX_PTS;
    V2 tmp[MAX_PTS];
    V2 chain[2 * MAX_PTS];
    memcpy(tmp, pts, (size_t)n * sizeof(V2));
    qsort(tmp, (size_t)n, sizeof(V2), cmp_v2);
    int m = 0;
    for (int i = 0; i < n; i++) {
        while (m >= 2 && cross2(chain[m - 2], chain[m - 1], tmp[i]) <= 0.0f) m--;
        chain[m++] = tmp[i];
    }
    int lower = m;
    for (int i = n - 2; i >= 0; i--) {
        while (m > lower && cross2(chain[m - 2], chain[m - 1], tmp[i]) <= 0.0f) m--;
        chain[m++] = tmp[i];
    }
    m--; /* closing point duplicates the first */
    if (m > n) m = n;
    if (m > MAX_PTS) m = MAX_PTS;
    memcpy(out, chain, (size_t)m * sizeof(V2));
    return m;
}

static float poly_area(const V2* p, int n) {
    if (n < 3) return 0.0f;
    float a = 0.0f;
    for (int i = 0; i < n; i++) {
        V2 u = p[i], v = p[(i + 1) % n];
        a += u.x * v.y - u.y * v.x;
    }
    return fabsf(a) * 0.5f;
}

static int point_in_poly(V2 p, const V2* poly, int n) {
    int inside = 0;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        V2 a = poly[j], b = poly[i];
        if (((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x))
            inside = !inside;
    }
    return inside;
}

/* ray vs convex polygon edges (any edge crossing at t >= 0) */
static int ray_poly(V2 o, V2 d, const V2* p, int n) {
    for (int i = 0; i < n; i++) {
        V2 a = p[i], b = p[(i + 1) % n];
        float ex = b.x - a.x, ey = b.y - a.y;
        float den = d.x * ey - d.y * ex;
        if (fabsf(den) < 1e-12f) continue;
        float qx = a.x - o.x, qy = a.y - o.y;
        float t = (qx * ey - qy * ex) / den;
        float u = (qx * d.y - qy * d.x) / den;
        if (t >= 0.0f && u >= 0.0f && u <= 1.0f) return 1;
    }
    return 0;
}

/* ---------------- criteria ----------------
 * sign = +1 correct : low concavity passes (C <= thr)
 * sign = -1 inverted: high "concavity" passes (C >= thr)
 */
static int concavity_gate(float C, float thr, int sign) {
    return sign > 0 ? (C <= thr) : (C >= thr);
}
static int tightness_gate(float T, float thr) { return T <= thr; }
static int hull_decision(float C, float T, int sign) {
    return concavity_gate(C, CONCAVE_THR, sign) &&
           tightness_gate(T, TIGHT_THR);
}

/* ---------------- clustering + proxy metrics ---------------- */
typedef struct {
    int idx[MAX_MEMBERS];
    int n;
    V2 mn, mx;
    V2 hull[MAX_PTS];
    int nh;
    float C;       /* actual concavity: 1 - covered / hull-area samples */
    float T;       /* mean nearest-neighbour centre dist / mean size */
    float weight;  /* hull area, for tier-level averaging */
    int use_hull_correct;
    int use_hull_inverted;
} Proxy;

static void cluster_metrics(const Obj* objs, const int* idx, int n, Proxy* p) {
    p->n = n;
    p->mn = (V2){ 1e30f, 1e30f };
    p->mx = (V2){ -1e30f, -1e30f };
    V2 pts[MAX_PTS];
    int np = 0;
    for (int k = 0; k < n; k++) {
        const Obj* o = &objs[idx[k]];
        for (int c = 0; c < 4; c++) {
            V2 q = o->corner[c];
            if (q.x < p->mn.x) p->mn.x = q.x;
            if (q.y < p->mn.y) p->mn.y = q.y;
            if (q.x > p->mx.x) p->mx.x = q.x;
            if (q.y > p->mx.y) p->mx.y = q.y;
            if (np < MAX_PTS) pts[np++] = q;
        }
    }
    p->nh = hull2d(pts, np, p->hull);
    float ha = poly_area(p->hull, p->nh);
    float ba = (p->mx.x - p->mn.x) * (p->mx.y - p->mn.y);
    p->weight = ha > 1e-6f ? ha : (ba > 1e-6f ? ba : 1e-6f);

    /* concavity: regular grid over the AABB, reject samples outside hull */
    int accepted = 0, covered = 0;
    float sx = p->mx.x - p->mn.x, sy = p->mx.y - p->mn.y;
    for (int ix = 0; ix < SAMPLE_GRID; ix++) {
        for (int iy = 0; iy < SAMPLE_GRID; iy++) {
            V2 q = { p->mn.x + (ix + 0.5f) / SAMPLE_GRID * sx,
                     p->mn.y + (iy + 0.5f) / SAMPLE_GRID * sy };
            if (p->nh >= 3 && !point_in_poly(q, p->hull, p->nh)) continue;
            accepted++;
            for (int k = 0; k < n; k++) {
                if (point_in_obb(q, &objs[idx[k]])) { covered++; break; }
            }
        }
    }
    p->C = accepted > 0 ? 1.0f - (float)covered / (float)accepted : 1.0f;

    /* tightness: mean NN centre distance divided by mean full size */
    if (n < 2) {
        p->T = 0.0f;
    } else {
        double nn_sum = 0.0, scale_sum = 0.0;
        int nn_cnt = 0;
        for (int a = 0; a < n; a++) {
            float best = 1e30f;
            for (int b = 0; b < n; b++) {
                if (a == b) continue;
                V2 d = vsub(objs[idx[a]].c, objs[idx[b]].c);
                float dd = sqrtf(d.x * d.x + d.y * d.y);
                if (dd < best) best = dd;
            }
            if (best < 1e29f) { nn_sum += best; nn_cnt++; }
            scale_sum += 2.0f * fmaxf(objs[idx[a]].hx, objs[idx[a]].hy);
        }
        float nnd = (float)(nn_sum / fmax(nn_cnt, 1));
        float scale = (float)(scale_sum / n);
        p->T = scale > 1e-6f ? nnd / scale : 0.0f;
    }

    /* clusters below 3 members gain nothing from a hull: always box */
    int eligible = (n >= 3 && p->nh >= 3);
    p->use_hull_correct  = eligible && hull_decision(p->C, p->T, +1);
    p->use_hull_inverted = eligible && hull_decision(p->C, p->T, -1);
}

static int build_proxies(const Obj* objs, int n, Proxy* ps, int maxp) {
    int parent[MAX_OBJ];
    for (int i = 0; i < n; i++) parent[i] = i;
    /* union-find clustering by centre distance */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            V2 d = vsub(objs[i].c, objs[j].c);
            if (d.x * d.x + d.y * d.y > CLUSTER_R * CLUSTER_R) continue;
            int a = i, b = j;
            while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
            while (parent[b] != b) { parent[b] = parent[parent[b]]; b = parent[b]; }
            if (a != b) parent[b] = a;
        }
    }
    int root_map[MAX_OBJ];
    for (int i = 0; i < MAX_OBJ; i++) root_map[i] = -1;
    int np = 0;
    for (int i = 0; i < n; i++) {
        int r = i;
        while (parent[r] != r) { parent[r] = parent[parent[r]]; r = parent[r]; }
        if (root_map[r] < 0) {
            if (np >= maxp) break;
            root_map[r] = np++;
            ps[root_map[r]].n = 0;
        }
        Proxy* p = &ps[root_map[r]];
        if (p->n < MAX_MEMBERS) p->idx[p->n++] = i;
    }
    for (int i = 0; i < np; i++) cluster_metrics(objs, ps[i].idx, ps[i].n, &ps[i]);
    return np;
}

/* ---------------- deterministic scene ---------------- */
#define ARENA_X 60.0f
#define ARENA_Z 40.0f

static const V2 PILES[9] = {
    { -34, -20 }, { 0, -20 }, { 34, -20 },
    { -34,   0 }, { 0,   0 }, { 34,   0 },
    { -34,  20 }, { 0,  20 }, { 34,  20 }
};

static int gen_scene(Obj* objs, float spacing, uint32_t seed) {
    Rng rng = { seed };
    int n = 0;
    for (int pi = 0; pi < 9; pi++) {
        float px = PILES[pi].x, pz = PILES[pi].y;
        float yaw = rng_range(&rng, 0.0f, 3.14159265f);
        float ca = cosf(yaw), sa = sinf(yaw);
        float a = rng_range(&rng, 9.5f, 11.5f);
        float b = rng_range(&rng, 2.0f, 2.8f);
        int nx = (int)(a / spacing) + 2;
        int nz = (int)(b / spacing) + 2;
        for (int ix = -nx; ix <= nx; ix++) {
            for (int iz = -nz; iz <= nz; iz++) {
                float lx = ix * spacing + rng_range(&rng, -0.16f, 0.16f) * spacing;
                float lz = iz * spacing + rng_range(&rng, -0.16f, 0.16f) * spacing;
                if ((lx / a) * (lx / a) + (lz / b) * (lz / b) > 1.0f) continue;
                if (n >= MAX_OBJ) return n;
                obj_init(&objs[n++],
                         px + lx * ca - lz * sa,
                         pz + lx * sa + lz * ca,
                         rng_range(&rng, 0.0f, 3.14159265f),
                         rng_range(&rng, 0.42f, 0.62f),
                         rng_range(&rng, 0.42f, 0.62f));
            }
        }
    }
    for (int s = 0; s < 10; s++) {
        if (n >= MAX_OBJ) break;
        obj_init(&objs[n++],
                 rng_range(&rng, -58.0f, 58.0f),
                 rng_range(&rng, -38.0f, 38.0f),
                 rng_range(&rng, 0.0f, 3.14159265f),
                 rng_range(&rng, 0.38f, 0.58f),
                 rng_range(&rng, 0.38f, 0.58f));
    }
    return n;
}

/* ---------------- rays + oracle ---------------- */
typedef struct { V2 o, d; } Ray;

static void gen_rays(uint32_t seed, Ray* rays) {
    Rng rng = { seed };
    for (int i = 0; i < BENCH_RAYS; i++) {
        int side = (int)rng_range(&rng, 0.0f, 4.0f);
        V2 o, d;
        if (side == 0) {
            o = (V2){ -ARENA_X - 2.0f, rng_range(&rng, -ARENA_Z, ARENA_Z) };
            d = (V2){ 1.0f, rng_range(&rng, -0.25f, 0.25f) };
        } else if (side == 1) {
            o = (V2){ ARENA_X + 2.0f, rng_range(&rng, -ARENA_Z, ARENA_Z) };
            d = (V2){ -1.0f, rng_range(&rng, -0.25f, 0.25f) };
        } else if (side == 2) {
            o = (V2){ rng_range(&rng, -ARENA_X, ARENA_X), -ARENA_Z - 2.0f };
            d = (V2){ rng_range(&rng, -0.25f, 0.25f), 1.0f };
        } else {
            o = (V2){ rng_range(&rng, -ARENA_X, ARENA_X), ARENA_Z + 2.0f };
            d = (V2){ rng_range(&rng, -0.25f, 0.25f), -1.0f };
        }
        float l = sqrtf(d.x * d.x + d.y * d.y);
        rays[i] = (Ray){ o, (V2){ d.x / l, d.y / l } };
    }
}

static int oracle_ray(const Ray* r, const Obj* objs, int nobj) {
    for (int i = 0; i < nobj; i++)
        if (ray_obb(r->o, r->d, &objs[i])) return 1;
    return 0;
}

static int proxy_ray_hit(const Ray* r, const Proxy* p, int inverted) {
    if (!ray_aabb(r->o, r->d, p->mn, p->mx)) return 0;
    if (!(inverted ? p->use_hull_inverted : p->use_hull_correct)) return 1;
    return ray_poly(r->o, r->d, p->hull, p->nh);
}

/* Unified frame core. Boolean collision queries: the first confirmed OBB
 * hit resolves the probe. `budget` caps narrow-phase OBB tests; pass a huge
 * budget for the uncapped demand / timing path. When the budget is gone,
 * the current and later probes are unserved; oracle-confirmed real hits
 * among them are returned as misses. */
typedef struct { long used; int misses; } FrameResult;

static FrameResult run_frame(const Ray* rays, const Obj* objs, int nobj,
                             const Proxy* ps, int np, int inverted, long budget) {
    FrameResult fr = { 0, 0 };
    for (int ri = 0; ri < BENCH_RAYS; ri++) {
        const Ray* r = &rays[ri];
        int resolved = 0;
        for (int pi = 0; pi < np && !resolved; pi++) {
            const Proxy* p = &ps[pi];
            if (!proxy_ray_hit(r, p, inverted)) continue;
            for (int k = 0; k < p->n; k++) {
                if (fr.used >= budget) break;
                fr.used++;
                if (ray_obb(r->o, r->d, &objs[p->idx[k]])) { resolved = 1; break; }
            }
        }
        if (!resolved && fr.used >= budget) {
            for (int rj = ri; rj < BENCH_RAYS; rj++)
                fr.misses += oracle_ray(&rays[rj], objs, nobj);
            break;
        }
    }
    return fr;
}

static long naive_frame_work(const Ray* rays, const Obj* objs, int nobj, long* sink) {
    long tests = 0, hits = 0;
    for (int ri = 0; ri < BENCH_RAYS; ri++) {
        for (int i = 0; i < nobj; i++) {
            tests++;
            if (ray_obb(rays[ri].o, rays[ri].d, &objs[i])) { hits++; break; }
        }
    }
    *sink += tests + hits;
    return tests;
}

/* Geometric false-negative check with an unlimited narrow phase. */
static long geometric_false_negatives(const Obj* objs, int nobj,
                                      const Proxy* ps, int np, int inverted) {
    long fn = 0;
    for (int f = 0; f < GEO_FRAMES; f++) {
        Ray rays[BENCH_RAYS];
        gen_rays(9000u + (uint32_t)f * 131u, rays);
        for (int ri = 0; ri < BENCH_RAYS; ri++) {
            if (!oracle_ray(&rays[ri], objs, nobj)) continue;
            int hit = 0;
            for (int pi = 0; pi < np; pi++)
                if (proxy_ray_hit(&rays[ri], &ps[pi], inverted)) { hit = 1; break; }
            if (!hit) fn++;
        }
    }
    return fn;
}

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static double time_proxy_us(const Obj* objs, int nobj, const Proxy* ps, int np,
                           int inverted, long* sink) {
    for (int f = 0; f < 30; f++) {
        Ray rays[BENCH_RAYS];
        gen_rays(50000u + (uint32_t)f * 977u, rays);
        FrameResult fr = run_frame(rays, objs, nobj, ps, np, inverted, (1L << 60));
        *sink += fr.used;
    }
    double best = 1e30;
    for (int rep = 0; rep < TIME_REPS; rep++) {
        long local = 0;
        double t0 = now_sec();
        for (int f = 0; f < TIME_FRAMES; f++) {
            Ray rays[BENCH_RAYS];
            gen_rays((uint32_t)(rep * 100000 + f * 977 + 7), rays);
            FrameResult fr = run_frame(rays, objs, nobj, ps, np, inverted, (1L << 60));
            local += fr.used;
        }
        double dt = now_sec() - t0;
        if (dt < best) best = dt;
        *sink += local;
    }
    return best / TIME_FRAMES * 1e6;
}

static double time_naive_us(const Obj* objs, int nobj, long* sink) {
    for (int f = 0; f < 10; f++) {
        Ray rays[BENCH_RAYS];
        gen_rays(50000u + (uint32_t)f * 977u, rays);
        naive_frame_work(rays, objs, nobj, sink);
    }
    double best = 1e30;
    for (int rep = 0; rep < TIME_REPS; rep++) {
        long local = 0;
        double t0 = now_sec();
        for (int f = 0; f < TIME_FRAMES; f++) {
            Ray rays[BENCH_RAYS];
            gen_rays((uint32_t)(rep * 100000 + f * 977 + 7), rays);
            naive_frame_work(rays, objs, nobj, &local);
        }
        double dt = now_sec() - t0;
        if (dt < best) best = dt;
        *sink += local;
    }
    return best / TIME_FRAMES * 1e6;
}

/* ---------------- experiment row ---------------- */
typedef struct {
    const char* label;
    int nobj, np;
    int hull_c, box_c, hull_i, box_i;
    int tier_gate_c, tier_gate_i;
    float C, T;
    long miss_c, miss_i;
    int worst_c, worst_i;
    long narrow_c0, narrow_i0;
    long fn_c, fn_i;
    double us_c, us_i, us_naive;
} Tier;

static void run_tier(Tier* t, float spacing, uint32_t seed, long* sink) {
    static Obj objs[MAX_OBJ];
    static Proxy ps[MAX_PROXY];
    memset(ps, 0, sizeof(ps));
    t->nobj = gen_scene(objs, spacing, seed);
    t->np = build_proxies(objs, t->nobj, ps, MAX_PROXY);

    double cnum = 0, wsum = 0, tsum = 0;
    int tcnt = 0;
    t->hull_c = t->hull_i = 0;
    for (int i = 0; i < t->np; i++) {
        cnum += ps[i].C * ps[i].weight;
        wsum += ps[i].weight;
        if (ps[i].n >= 2) { tsum += ps[i].T; tcnt++; }
        if (ps[i].use_hull_correct) t->hull_c++;
        if (ps[i].use_hull_inverted) t->hull_i++;
    }
    t->C = (float)(cnum / wsum);
    t->T = tcnt ? (float)(tsum / tcnt) : 0.0f;
    t->box_c = t->np - t->hull_c;
    t->box_i = t->np - t->hull_i;
    t->tier_gate_c = hull_decision(t->C, t->T, +1);
    t->tier_gate_i = hull_decision(t->C, t->T, -1);

    t->miss_c = t->miss_i = 0;
    t->worst_c = t->worst_i = 0;
    for (int f = 0; f < BENCH_FRAMES; f++) {
        Ray rays[BENCH_RAYS];
        gen_rays(9000u + (uint32_t)f * 131u + seed, rays);
        FrameResult fc = run_frame(rays, objs, t->nobj, ps, t->np, 0, NARROW_BUDGET);
        FrameResult fi = run_frame(rays, objs, t->nobj, ps, t->np, 1, NARROW_BUDGET);
        t->miss_c += fc.misses;
        t->miss_i += fi.misses;
        if (fc.misses > t->worst_c) t->worst_c = fc.misses;
        if (fi.misses > t->worst_i) t->worst_i = fi.misses;
    }
    t->fn_c = geometric_false_negatives(objs, t->nobj, ps, t->np, 0);
    t->fn_i = geometric_false_negatives(objs, t->nobj, ps, t->np, 1);
    {
        Ray rays0[BENCH_RAYS];
        gen_rays(9000u + seed, rays0); /* benchmark frame 0 */
        t->narrow_c0 = run_frame(rays0, objs, t->nobj, ps, t->np, 0, (1L << 60)).used;
        t->narrow_i0 = run_frame(rays0, objs, t->nobj, ps, t->np, 1, (1L << 60)).used;
    }
    t->us_c = time_proxy_us(objs, t->nobj, ps, t->np, 0, sink);
    t->us_i = time_proxy_us(objs, t->nobj, ps, t->np, 1, sink);
    t->us_naive = time_naive_us(objs, t->nobj, sink);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("check_%-34s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    printf("# proxy: distance-clustered collision proxies (AABB | convex hull)\n");
    printf("model=2D top-down OBBs, monotone-chain hull, %d ray probes/frame\n", BENCH_RAYS);
    printf("cluster_r=%.2f  criteria: C<=%.2f AND T<=%.2f -> hull else AABB\n",
           CLUSTER_R, CONCAVE_THR, TIGHT_THR);
    printf("narrow_budget=%ld tests/frame, miss frames=%d, timing best-of-%d x %d frames\n\n",
           NARROW_BUDGET, BENCH_FRAMES, TIME_REPS, TIME_FRAMES);

    Tier tiers[3];
    const float spacings[3] = { 3.2f, 1.75f, 1.35f };
    const uint32_t seeds[3] = { 101u, 202u, 303u };
    const char* labels[3] = { "low", "mid", "high" };
    long sink = 0;
    for (int i = 0; i < 3; i++) {
        tiers[i].label = labels[i];
        run_tier(&tiers[i], spacings[i], seeds[i], &sink);
    }

    /* ---- headline: three tiers x (proxies, us/frame, misses) = 9 numbers ---- */
    printf("== nine numbers, correct concavity criterion ==\n");
    printf("density objects proxies hull/box C_act T_act  proxy_us naive_us speedup misses/%df worst\n",
           BENCH_FRAMES);
    for (int i = 0; i < 3; i++) {
        Tier* t = &tiers[i];
        printf("%-5s %4d %7d %4d/%-3d %5.3f %5.2f %9.1f %8.1f %6.2fx %6ld %5d\n",
               t->label, t->nobj, t->np, t->hull_c, t->box_c, t->C, t->T,
               t->us_c, t->us_naive, t->us_naive / t->us_c, t->miss_c, t->worst_c);
    }
    printf("NINE_proxies=%d,%d,%d\n", tiers[0].np, tiers[1].np, tiers[2].np);
    printf("NINE_collision_us_per_frame=%.1f,%.1f,%.1f\n",
           tiers[0].us_c, tiers[1].us_c, tiers[2].us_c);
    printf("NINE_misses=%ld,%ld,%ld\n\n",
           tiers[0].miss_c, tiers[1].miss_c, tiers[2].miss_c);

    /* ---- requested side-by-side: actual concavity vs selected shape ---- */
    printf("== actual concavity/tightness vs selected proxy shape ==\n");
    printf("tier C_act T_act C_gate T_gate | correct tier-shape (hull/box) | inverted tier-shape (hull/box)\n");
    for (int i = 0; i < 3; i++) {
        Tier* t = &tiers[i];
        const char* sc = t->tier_gate_c ? "hull" : "bbox";
        const char* si = t->tier_gate_i ? "hull" : "bbox";
        printf("%-5s %.3f %.2f %-6s %-6s | %-4s (%2d/%-2d)            | %-4s (%2d/%-2d)\n",
               t->label, t->C, t->T,
               concavity_gate(t->C, CONCAVE_THR, +1) ? "pass" : "fail",
               tightness_gate(t->T, TIGHT_THR) ? "pass" : "fail",
               sc, t->hull_c, t->box_c, si, t->hull_i, t->box_i);
    }
    printf("\n");

    /* ---- inverted-sign control miss table ---- */
    printf("== inverted concavity control (C >= %.2f AND T <= %.2f -> hull) ==\n",
           CONCAVE_THR, TIGHT_THR);
    printf("density hull/box proxy_us misses/%df worst\n", BENCH_FRAMES);
    for (int i = 0; i < 3; i++) {
        Tier* t = &tiers[i];
        printf("%-5s %4d/%-3d %9.1f %6ld %5d\n",
               t->label, t->hull_i, t->box_i, t->us_i, t->miss_i, t->worst_i);
    }
    int inv_high_max = tiers[2].miss_i > tiers[1].miss_i &&
                       tiers[2].miss_i > tiers[0].miss_i && tiers[2].miss_i > 0;
    int corr_high_max = tiers[2].miss_c > tiers[1].miss_c ||
                        tiers[2].miss_c > tiers[0].miss_c;
    printf("diagnostic_inverted_high_density_is_worst=%d\n", inv_high_max);
    printf("diagnostic_correct_high_density_is_worst=%d\n\n", corr_high_max);

    printf("== narrow-phase demand on benchmark frame 0 (correct / inverted) ==\n");
    for (int i = 0; i < 3; i++)
        printf("%-5s narrow_tests=%ld / %ld\n",
               tiers[i].label, tiers[i].narrow_c0, tiers[i].narrow_i0);
    printf("\n");

    /* ---------------- self-checks ---------------- */
    V2 sq[4] = { {0,0},{1,0},{1,1},{0,1} };
    V2 hs[8];
    int nh = hull2d(sq, 4, hs);
    check("hull_square_4verts", nh == 4);
    check("hull_square_area_one", fabsf(poly_area(hs, nh) - 1.0f) < 1e-5f);

    Obj box; obj_init(&box, 0, 0, 0.0f, 1.0f, 1.0f);
    check("ray_hits_box", ray_obb((V2){ -3, 0.2f }, (V2){ 1, 0 }, &box) == 1);
    check("ray_misses_box", ray_obb((V2){ -3, 2.0f }, (V2){ 1, 0 }, &box) == 0);

    check("criteria_both_pass_hull", hull_decision(0.1f, 1.0f, +1) == 1);
    check("criteria_high_concavity_box", hull_decision(0.9f, 1.0f, +1) == 0);
    check("criteria_loose_arrangement_box", hull_decision(0.1f, 9.0f, +1) == 0);
    check("criteria_are_AND", hull_decision(0.9f, 9.0f, +1) == 0);
    check("inverted_flips_concavity_gate",
          hull_decision(0.9f, 1.0f, -1) == 1 && hull_decision(0.1f, 1.0f, -1) == 0);

    /* sign sanity on controlled clusters: filled block vs sparse ring */
    {
        Obj fo[64]; int fi[64]; int nf = 0;
        for (int ix = 0; ix < 6; ix++)
            for (int iz = 0; iz < 6; iz++) {
                obj_init(&fo[nf], (ix - 2.5f) * 0.9f, (iz - 2.5f) * 0.9f, 0.3f, 0.34f, 0.34f);
                fi[nf] = nf; nf++;
            }
        Proxy fp; memset(&fp, 0, sizeof fp);
        cluster_metrics(fo, fi, nf, &fp);

        Obj ro[16]; int ri[16]; int nr = 0;
        for (int k = 0; k < 10; k++) {
            float a = k * 6.2831853f / 10.0f;
            obj_init(&ro[nr], 5.0f * cosf(a), 5.0f * sinf(a), a, 0.34f, 0.34f);
            ri[nr] = nr; nr++;
        }
        Proxy rp; memset(&rp, 0, sizeof rp);
        cluster_metrics(ro, ri, nr, &rp);
        check("concavity_sparse_gt_filled", rp.C > fp.C + 0.25f);
        check("concavity_gate_direction_sane",
              concavity_gate(fp.C, CONCAVE_THR, +1) == 1 &&
              concavity_gate(rp.C, CONCAVE_THR, +1) == 0);
    }

    check("tier_concavity_monotone",
          tiers[0].C > tiers[1].C + 0.03f && tiers[1].C > tiers[2].C + 0.03f);
    check("tier_tightness_monotone", tiers[0].T > tiers[1].T && tiers[1].T > tiers[2].T);
    check("proxy_count_compresses_objects",
          tiers[0].np * 2 < tiers[0].nobj &&
          tiers[1].np * 2 < tiers[1].nobj &&
          tiers[2].np * 2 < tiers[2].nobj);
    check("high_density_correct_picks_hulls", tiers[2].hull_c >= 9);

    int geo_ok = 1, zero_miss_c = 1;
    for (int i = 0; i < 3; i++) {
        if (tiers[i].fn_c != 0 || tiers[i].fn_i != 0) geo_ok = 0;
        if (tiers[i].miss_c != 0) zero_miss_c = 0;
    }
    check("conservative_proxies_no_geometric_FN", geo_ok);
    check("correct_budget_meets_zero_misses", zero_miss_c);
    check("diag_inverted_high_density_worst", inv_high_max);
    check("diag_correct_rule_not_high_worst", !corr_high_max);
    check("inverted_box_wastes_narrow_phase",
          tiers[2].narrow_i0 > tiers[2].narrow_c0 * 1.10f);
    check("timing_measured", tiers[2].us_c > 1.0 && sink != 0);

    printf("\nall_ok=%d\n", failures == 0 ? 1 : 0);
    return failures == 0 ? 0 : 1;
}
