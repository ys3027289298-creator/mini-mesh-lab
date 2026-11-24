
/* simplify: edge-collapse simplification with convex-edge vertex locking.
 *
 * Collision mesh simplification rounds convex step edges and the walking
 * character falls through. A locked vertex blocks every collapse incident to
 * it. Which vertices get locked is decided by two coupled indicators:
 *   - dihedral angle at the edge (sharpness)
 *   - adjacent face size (a big riser/tread pair loses more support when cut)
 * Sweeping the lock budget gives the lock-ratio vs height-loss frontier and
 * the lock-ratio vs reachable reduction frontier; the target reduction is
 * inverted into a lock ratio, then three staircases of different step height
 * are simplified and the fall counts compared. An angle-only control run
 * shows what ignoring the size term costs on the tallest staircase.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NSTEP 10         /* steps per staircase */
#define SUBD  3          /* subdivisions per profile segment */
#define NSEG  (2 * NSTEP + 2)         /* tread+riser per step, tread of the
                                       * last step, plus a landing extension */
#define EMAX  (NSEG * SUBD)             /* profile edges            */
#define PMAX  (EMAX + 1)                /* profile vertices         */
#define MAXV  (2 * PMAX + 16)           /* extruded mesh vertices   */
#define MAXF  (2 * EMAX + 16)           /* extruded mesh triangles  */
#define SWEEPK (2 * NSTEP + 4)          /* sweep 0..24 locks        */

#define WIDTH   1.0f
#define DEPTH   1.0f                    /* tread depth d            */
#define FALL_DL 0.30f                   /* support loss counted as a fall */
#define R_TARGET 0.75f                  /* wanted simplification rate */

typedef struct { float x, y; } P2;
typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

/* tag values for profile vertices */
#define TAG_STRAIGHT 0
#define TAG_CONVEX   1
#define TAG_CONCAVE (-1)

typedef struct {
    float h[NSTEP];
    P2    prof[PMAX];
    int   tag[PMAX];
    int   np;
    V3    v[MAXV];
    int   nv;
    Tri   f[MAXF];
    int   nf;
} Stair;

/* ---------- build: 2D walking profile, then extrude along z ---------- */

static void build_stair(Stair* s, const float* heights) {
    memset(s, 0, sizeof(*s));
    memcpy(s->h, heights, sizeof(float) * NSTEP);

    int p = 0;
    s->prof[p] = (P2){0.0f, 0.0f};
    s->tag[p] = TAG_STRAIGHT;
    p++;
    float x = 0.0f;
    for (int i = 0; i < NSTEP; i++) {
        float y0 = (i == 0) ? 0.0f : heights[i - 1];
        float y1 = heights[i];
        /* riser segment: vertical, ends in a concave (inner) corner */
        for (int k = 1; k <= SUBD; k++) {
            float t = (float)k / SUBD;
            s->prof[p].x = x;
            s->prof[p].y = y0 + (y1 - y0) * t;
            s->tag[p] = (k == SUBD) ? TAG_CONCAVE : TAG_STRAIGHT;
            p++;
        }
        /* tread segment: horizontal, ends in a convex (outer) corner */
        for (int k = 1; k <= SUBD; k++) {
            float t = (float)k / SUBD;
            s->prof[p].x = x + DEPTH * t;
            s->prof[p].y = y1;
            s->tag[p] = (k == SUBD) ? TAG_CONVEX : TAG_STRAIGHT;
            p++;
        }
        x += DEPTH;
    }
    /* landing extension: two more tread-depth segments beyond the last step's
     * tread, so the final convex edge is an interior edge with a quad on each
     * side (its own tread behind, the landing ahead), not a mesh boundary */
    for (int seg = 0; seg < 2; seg++) {
        for (int k = 1; k <= SUBD; k++) {
            float t = (float)k / SUBD;
            s->prof[p].x = x + DEPTH * (seg + t);
            s->prof[p].y = heights[NSTEP - 1];
            s->tag[p] = TAG_STRAIGHT;
            p++;
        }
    }
    s->np = p;

    /* extrude: two z layers, vertex (profile i, layer l) = 2*i + l */
    for (int i = 0; i < s->np; i++) {
        s->v[2 * i]     = (V3){s->prof[i].x, s->prof[i].y, 0.0f};
        s->v[2 * i + 1] = (V3){s->prof[i].x, s->prof[i].y, WIDTH};
    }
    s->nv = 2 * s->np;
    for (int e = 0; e < s->np - 1; e++) {
        int a = 2 * e, b = 2 * (e + 1);
        s->f[s->nf++] = (Tri){a, b, b + 1};
        s->f[s->nf++] = (Tri){a, b + 1, a + 1};
    }
}

/* ---------- mesh measurements: face area, dihedral, adjacent size ---------- */

static V3 vsub(V3 a, V3 b) { return (V3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static V3 vcross(V3 a, V3 b) {
    return (V3){a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
static float vdot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static float vlen(V3 a) { return sqrtf(vdot(a, a)); }
static V3 vunit(V3 a) {
    float l = vlen(a);
    return l > 1e-12f ? (V3){a.x/l, a.y/l, a.z/l} : (V3){0,0,0};
}

/* incident full-segment quad area for each profile vertex. Each profile
 * vertex lies on two whole segments (riser/tread); measure each segment
 * from its end vertices (segment j spans j*SUBD .. (j+1)*SUBD) and add it
 * to every subdivision vertex on that segment. */
static void adjacent_sizes(const Stair* s, float* area) {
    for (int i = 0; i < s->np; i++) area[i] = 0.0f;
    for (int j = 0; j < NSEG; j++) {
        int e0 = j * SUBD, e1 = (j + 1) * SUBD;
        float dx = s->prof[e1].x - s->prof[e0].x;
        float dy = s->prof[e1].y - s->prof[e0].y;
        float qa = sqrtf(dx * dx + dy * dy) * WIDTH;
        for (int i = e0; i <= e1; i++) area[i] += qa;
    }
}

/* dihedral at profile vertex i: angle between the normals of the two
 * subdivision quads meeting there. The polyline edges touching vertex i are
 * edges i-1 and i; extrusion quad on edge e holds triangles 2*e and 2*e+1. */
static float dihedral_at(const Stair* s, int i) {
    if (i <= 0 || i >= s->np - 1) return 0.0f;
    V3 n0 = {0,0,0}, n1 = {0,0,0};
    int q0 = 2 * (i - 1), q1 = 2 * i;           /* first triangle of each quad */
    V3 t0a = vcross(vsub(s->v[s->f[q0].b], s->v[s->f[q0].a]),
                    vsub(s->v[s->f[q0].c], s->v[s->f[q0].a]));
    V3 t0b = vcross(vsub(s->v[s->f[q0+1].b], s->v[s->f[q0+1].a]),
                    vsub(s->v[s->f[q0+1].c], s->v[s->f[q0+1].a]));
    V3 t1a = vcross(vsub(s->v[s->f[q1].b], s->v[s->f[q1].a]),
                    vsub(s->v[s->f[q1].c], s->v[s->f[q1].a]));
    V3 t1b = vcross(vsub(s->v[s->f[q1+1].b], s->v[s->f[q1+1].a]),
                    vsub(s->v[s->f[q1+1].c], s->v[s->f[q1+1].a]));
    n0 = vunit((V3){t0a.x+t0b.x, t0a.y+t0b.y, t0a.z+t0b.z});
    n1 = vunit((V3){t1a.x+t1b.x, t1a.y+t1b.y, t1a.z+t1b.z});
    float c = vdot(n0, n1);
    if (c > 1.0f) c = 1.0f;
    if (c < -1.0f) c = -1.0f;
    return acosf(c);
}

/* ---------- lock selection ---------- */

/* use_size: gate-then-rank selection. The dihedral angle is the gate: a
 * sharper edge is never skipped in favour of a blunter one, so the budget
 * works through angle bands from sharpest down. Within one band the adjacent
 * face size is the sole ranking key (a big riser/tread junction loses more
 * support when cut, so it takes the lock first). An angle-only rule instead
 * breaks in-band ties by vertex index with no size information. This matters
 * for a convex edge that is coplanar with its landing (dihedral 0): under
 * angle-only it is unreachable behind every 90-degree corner, while the size
 * band lets it be locked once the sharper bands are spent. */
static void pick_locks(const Stair* s, int k, int use_size, int* lock) {
    float area[PMAX];
    adjacent_sizes(s, area);

    for (int i = 0; i < s->np; i++) lock[i] = 0;
    for (int n = 0; n < k; n++) {
        int best = -1;
        float best_angle = -1e30f, best_size = -1e30f;
        for (int i = 0; i < s->np; i++) {
            if (lock[i] || s->tag[i] == TAG_STRAIGHT) continue;
            float ang = dihedral_at(s, i) / (float)M_PI;
            int better;
            if (use_size) {
                /* sharpest non-empty angle band first, then largest face */
                better = (ang > best_angle + 1e-9f) ||
                         (fabsf(ang - best_angle) <= 1e-9f &&
                          (area[i] > best_size + 1e-9f ||
                           (fabsf(area[i] - best_size) <= 1e-9f &&
                            (best < 0 || i < best))));
            } else {
                /* angle only: sharper wins, index breaks in-band ties */
                better = (ang > best_angle + 1e-9f) ||
                         (fabsf(ang - best_angle) <= 1e-9f &&
                          (best < 0 || i < best));
            }
            if (better) {
                best = i; best_angle = ang;
                if (use_size) best_size = area[i];
            }
        }
        if (best < 0) break;              /* no more corner candidates */
        lock[best] = 1;
    }
}

/* ---------- greedy edge collapse on the profile polyline ---------- */

static int simplify_profile(Stair* s, const int* lock, P2* out, int* alive) {
    for (int i = 0; i < s->np; i++) { out[i] = s->prof[i]; alive[i] = 1; }
    int collapses = 0;
    for (;;) {
        int best = -1;
        float best_cost = 1e30f;
        for (int e = 1; e < s->np - 1; e++) {
            if (!alive[e] || lock[e]) continue;
            int lo = e - 1; while (!alive[lo]) lo--;
            int hi = e + 1; while (!alive[hi]) hi++;
            if (lock[lo] && lock[hi]) continue;      /* locked ends must not move */
            /* polyline simplification cost: perpendicular deviation of e from
             * the chord joining its alive neighbours (zero on straight runs) */
            float ax = out[hi].x - out[lo].x;
            float ay = out[hi].y - out[lo].y;
            float alen2 = ax * ax + ay * ay;
            float px = out[e].x - out[lo].x;
            float py = out[e].y - out[lo].y;
            float cross2 = ax * py - ay * px;
            float cost = (cross2 * cross2) / (alen2 > 1e-18f ? alen2 : 1e-18f);
            if (cost < best_cost) { best_cost = cost; best = e; }
        }
        if (best < 0) break;
        int lo = best - 1; while (!alive[lo]) lo--;
        int hi = best + 1; while (!alive[hi]) hi++;
        if (lock[lo]) out[best] = out[hi];          /* keep locked endpoint fixed */
        else if (lock[hi]) out[best] = out[lo];
        else {
            out[best].x = 0.5f * (out[lo].x + out[hi].x);
            out[best].y = 0.5f * (out[lo].y + out[hi].y);
        }
        alive[best] = 0;
        collapses++;
    }
    return collapses;
}

/* support height of the simplified surface at query x: the character stands
 * on top, so take the upper envelope — the maximum y over every surviving
 * segment that spans qx (a collapsed concave corner leaves a diagonal where
 * the tread used to be; the tread edge above it is what still supports the
 * foot). The alive flags are sparse, so first compact the surviving points
 * into a contiguous polyline; adjacent indices in that array are the real
 * neighbours after collapses. */
static float support_y(const Stair* s, const P2* pos, const int* alive, float qx) {
    P2 chain[PMAX];
    int nc = 0;
    for (int i = 0; i < s->np; i++)
        if (alive[i]) chain[nc++] = pos[i];

    float hi_y = 0.0f;
    int found = 0;
    for (int i = 0; i < nc - 1; i++) {
        if (qx < chain[i].x - 1e-9f || qx > chain[i + 1].x + 1e-9f) continue;
        float dx = chain[i + 1].x - chain[i].x;
        if (dx <= 1e-12f) continue;                 /* vertical riser: no floor */
        float t = (qx - chain[i].x) / dx;
        float y = chain[i].y + t * (chain[i + 1].y - chain[i].y);
        if (!found || y > hi_y) { hi_y = y; found = 1; }
    }
    return hi_y;
}

/* probe every convex corner at the middle of its own tread. The convex
 * corner ends the tread; the tread begins at the preceding concave corner,
 * one profile segment back. Querying the tread midpoint (not a fixed offset
 * from the edge) guarantees the probe sits over tread geometry rather than
 * the riser, so a rounded edge that slopes the tread is what gets measured. */
static void probe_losses(const Stair* s, const P2* pos, const int* alive,
                         float* loss, int* nfall) {
    *nfall = 0;
    for (int i = 0; i < s->np; i++) loss[i] = 0.0f;
    for (int i = 0; i < s->np; i++) {
        if (s->tag[i] != TAG_CONVEX) continue;
        int tread_start = i - SUBD;                 /* concave corner before it */
        float qx = 0.5f * (s->prof[tread_start].x + s->prof[i].x);
        float drop = s->prof[i].y - support_y(s, pos, alive, qx);
        if (drop < 0.0f) drop = 0.0f;
        loss[i] = drop;
        if (drop > FALL_DL) (*nfall)++;
    }
}

static float max_loss(const Stair* s, const float* loss) {
    float m = 0.0f;
    for (int i = 0; i < s->np; i++)
        if (s->tag[i] == TAG_CONVEX && loss[i] > m) m = loss[i];
    return m;
}

/* ---------- one full simplification experiment ---------- */

typedef struct {
    int   locks_profile;      /* locked profile vertices */
    int   locks_mesh;         /* locked 3D vertices (2 layers) */
    float lock_ratio;
    float max_height_loss;
    float reduction;
    int   falls;
} Result;

static Result run_once(Stair* s, int k, int use_size) {
    int lock[PMAX];
    pick_locks(s, k, use_size, lock);
    P2 pos[PMAX];
    int alive[PMAX];
    int collapses = simplify_profile(s, lock, pos, alive);
    float loss[PMAX];
    int falls;
    probe_losses(s, pos, alive, loss, &falls);
    Result r;
    r.locks_profile = 0;
    for (int i = 0; i < s->np; i++) r.locks_profile += lock[i];
    r.locks_mesh = 2 * r.locks_profile;
    r.lock_ratio = (float)r.locks_profile / s->np;
    r.max_height_loss = max_loss(s, loss);
    r.reduction = (float)collapses / EMAX;
    r.falls = falls;
    return r;
}

/* ---------- test harness ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    static const float hS[NSTEP] = {0.15f,0.30f,0.45f,0.60f,0.75f,0.90f,1.05f,1.20f,1.35f,1.50f};
    static const float hM[NSTEP] = {0.30f,0.45f,0.60f,0.75f,0.90f,1.05f,1.20f,1.35f,1.50f,1.65f};
    /* tall staircase with increasing riser heights so the size term spreads:
     * convex-corner adjacent area = 1.0 (tread) + this riser, 0.30 -> 0.75 */
    static const float hT[NSTEP] = {0.30f,0.75f,1.35f,2.10f,3.00f,4.05f,5.25f,6.60f,8.10f,9.75f};
    const float* hh[3] = { hS, hM, hT };
    const char* names[3] = { "S_low", "M_mid", "T_tall" };

    Stair meshes[3];
    for (int g = 0; g < 3; g++) build_stair(&meshes[g], hh[g]);

    printf("## params\n");
    printf("steps=%d subdiv=%d profile_edges=%d profile_verts=%d mesh_verts=%d mesh_tris=%d\n",
           NSTEP, SUBD, EMAX, PMAX, meshes[0].nv, meshes[0].nf);
    printf("width=%.2f depth=%.2f scoring=gate_angle_rank_size fall_delta=%.2f target_reduction=%.2f\n",
           WIDTH, DEPTH, FALL_DL, R_TARGET);

    /* ---- unit checks ---- */
    printf("## unit_checks\n");
    Stair* t = &meshes[2];
    int nconvex = 0, nconcave = 0;
    for (int i = 0; i < t->np; i++) {
        if (t->tag[i] == TAG_CONVEX) nconvex++;
        if (t->tag[i] == TAG_CONCAVE) nconcave++;
    }
    check("build_counts", t->np == PMAX && t->nv == 2 * PMAX && t->nf == 2 * EMAX);
    check("corner_counts", nconvex == NSTEP && nconcave == NSTEP);

    int first_convex = -1;
    for (int i = 0; i < t->np; i++) if (t->tag[i] == TAG_CONVEX) { first_convex = i; break; }
    check("convex_dihedral_90", fabsf(dihedral_at(t, first_convex) - (float)M_PI / 2) < 1e-4f);
    check("straight_dihedral_0", fabsf(dihedral_at(t, first_convex - 1)) < 1e-4f);

    float area[PMAX];
    adjacent_sizes(t, area);
    /* tallest corner: adjacent quads are riser (h * W) and tread (d * W) */
    int last_convex = -1;
    for (int i = 0; i < t->np; i++) if (t->tag[i] == TAG_CONVEX) last_convex = i;
    /* last convex corner sits between the final tread and the landing */
    check("last_corner_two_treads", fabsf(area[last_convex] - 2.0f * DEPTH * WIDTH) < 1e-4f);
    /* an interior convex corner sits between its riser and its tread */
    int mid_convex = (2 * (NSTEP / 2) + 1) * SUBD;
    float mid_riser = t->h[NSTEP / 2] - t->h[NSTEP / 2 - 1];
    check("adjacent_area_size_link",
          fabsf(area[mid_convex] - (mid_riser + DEPTH * WIDTH)) < 1e-4f);
    check("area_grows_with_step", area[mid_convex] - area[first_convex] > 0.25f);

    /* every lock must land on a step-edge corner, never a straight point */
    int probe_lock[PMAX];
    pick_locks(t, NSTEP / 2, 1, probe_lock);
    int only_corner = 1, nlocked = 0;
    for (int i = 0; i < t->np; i++) {
        if (!probe_lock[i]) continue;
        nlocked++;
        if (t->tag[i] == TAG_STRAIGHT) only_corner = 0;
    }
    check("locks_only_on_corners", only_corner && nlocked == NSTEP / 2);

    /* no locks: every interior vertex collapses, only the two polyline
     * endpoints remain, so EMAX-1 of EMAX edges are removable */
    Result none = run_once(t, 0, 1);
    check("no_lock_full_reduction",
          fabsf(none.reduction - (float)(EMAX - 1) / EMAX) < 1e-6f);
    check("no_lock_causes_falls", none.falls >= 1 && none.max_height_loss > FALL_DL);

    /* all step-edge corners locked: treads and risers stay intact, zero loss */
    Result all = run_once(t, 2 * NSTEP, 1);
    check("all_corners_locked_zero_loss", all.max_height_loss == 0.0f && all.falls == 0);
    check("locks_cap_reduction", all.reduction < 1.0f);

    /* ---- frontier sweep: L(p) and R_max(p) ---- */
    printf("## frontier_sweep combined_scoring\n");
    printf("k locks_profile locks_3d lock_ratio  L_S      L_M      L_T      Rmax_S   Rmax_M   Rmax_T\n");
    float rmax[SWEEPK + 1][3];
    for (int k = 0; k <= SWEEPK; k++) {
        Result r[3];
        for (int g = 0; g < 3; g++) r[g] = run_once(&meshes[g], k, 1);
        for (int g = 0; g < 3; g++) rmax[k][g] = r[g].reduction;
        printf("%-2d %-13d %-9d %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f\n",
               k, r[0].locks_profile, r[0].locks_mesh, r[0].lock_ratio,
               r[0].max_height_loss, r[1].max_height_loss, r[2].max_height_loss,
               r[0].reduction, r[1].reduction, r[2].reduction);
    }

    /* inversion: largest lock budget whose reachable reduction still meets target */
    int kstar = -1;
    for (int k = 0; k <= SWEEPK; k++) {
        int ok_all = 1;
        for (int g = 0; g < 3; g++) if (rmax[k][g] + 1e-9f < R_TARGET) ok_all = 0;
        if (ok_all) kstar = k;
    }
    float pstar = (float)kstar / PMAX;
    printf("## inversion target_reduction=%.2f\n", R_TARGET);
    printf("k_star=%d p_star=%.4f locks_3d=%d\n", kstar, pstar, 2 * kstar);
    check("inversion_found", kstar > 0 && kstar < NSTEP);

    int mono = 1;
    for (int k = 1; k <= SWEEPK; k++)
        for (int g = 0; g < 3; g++)
            if (rmax[k][g] > rmax[k - 1][g] + 1e-9f) mono = 0;
    check("rmax_nonincreasing", mono);
    check("inversion_boundary",
          rmax[kstar][0] + 1e-9f >= R_TARGET && rmax[kstar + 1][0] < R_TARGET);

    /* ---- main runs: combined scoring at the inverted lock ratio ---- */
    printf("## main_runs combined angle+size at p_star\n");
    printf("group     locks_3d lock_ratio height_loss reduction falls\n");
    Result mainr[3];
    for (int g = 0; g < 3; g++) {
        mainr[g] = run_once(&meshes[g], kstar, 1);
        printf("%-9s %-9d %.4f     %.4f      %.4f    %d\n",
               names[g], mainr[g].locks_mesh, mainr[g].lock_ratio,
               mainr[g].max_height_loss, mainr[g].reduction, mainr[g].falls);
    }
    /* combined scoring keeps the two lower staircases fall-free at the
     * inverted budget; the tall-one comparison lives after the control run */
    check("combined_low_mid_zero_falls", mainr[0].falls == 0 && mainr[1].falls == 0);

    /* ---- control: angle-only locks, same budget ---- */
    printf("## control_runs angle_only same_budget\n");
    printf("group     locks_3d lock_ratio height_loss reduction falls\n");
    Result ctlr[3];
    for (int g = 0; g < 3; g++) {
        ctlr[g] = run_once(&meshes[g], kstar, 0);
        printf("%-9s %-9d %.4f     %.4f      %.4f    %d\n",
               names[g], ctlr[g].locks_mesh, ctlr[g].lock_ratio,
               ctlr[g].max_height_loss, ctlr[g].reduction, ctlr[g].falls);
    }
    check("angle_only_tall_falls", ctlr[2].falls >= 1);
    check("combined_beats_angle_on_tall", mainr[2].falls < ctlr[2].falls);
    check("angle_only_worse_or_equal",
          ctlr[0].falls >= mainr[0].falls && ctlr[1].falls >= mainr[1].falls &&
          ctlr[2].falls >= mainr[2].falls);

    /* ---- diagnostic on the tallest staircase: angle vs size side by side ---- */
    printf("## diagnostic T_tall convex_corners dihedral_vs_adjacent_size\n");
    int lock_size[PMAX], lock_ang[PMAX];
    pick_locks(t, kstar, 1, lock_size);
    pick_locks(t, kstar, 0, lock_ang);
    P2 pos[PMAX]; int alive[PMAX]; float loss[PMAX]; int dummy_falls;
    simplify_profile(t, lock_ang, pos, alive);
    probe_losses(t, pos, alive, loss, &dummy_falls);
    printf("corner step_h dihedral_rad dihedral_deg adjacent_area lock_combined lock_angle_only support_loss fall\n");
    int ci = 0;
    for (int i = 0; i < t->np; i++) {
        if (t->tag[i] != TAG_CONVEX) continue;
        float ang = dihedral_at(t, i);
        printf("%-6d %.2f   %.4f       %.2f          %.4f         %d             %d              %.4f        %d\n",
               ci + 1, t->h[ci], ang, ang * 180.0f / (float)M_PI, area[i],
               lock_size[i], lock_ang[i], loss[i], loss[i] > FALL_DL);
        ci++;
    }
    check("size_spreads_risk", area[mid_convex] - area[first_convex] > 0.25f);
    int probe_size[PMAX], probe_ang[PMAX];
    pick_locks(t, NSTEP, 1, probe_size);
    pick_locks(t, NSTEP, 0, probe_ang);
    /* the last convex edge is coplanar with its landing, so its dihedral is
     * zero: no angle-based rule, combined or angle-only, can see it through
     * the angle gate. This is the precise failure mode of "lock by dihedral
     * alone". What distinguishes the edges that angle CAN see is face size:
     * within the 90-degree convex band the size ranking locks the biggest
     * riser junctions first, while angle-only spends the same budget on the
     * lowest-index edges. Demonstrate both facts. */
    check("flat_edge_invisible_to_angle",
          dihedral_at(t, last_convex) < 1e-4f &&
          probe_size[last_convex] == 0 && probe_ang[last_convex] == 0);
    int band_size[PMAX], band_ang[PMAX];
    pick_locks(t, NSTEP / 2, 1, band_size);
    pick_locks(t, NSTEP / 2, 0, band_ang);
    float locked_size_sum = 0.0f, locked_ang_sum = 0.0f;
    for (int i = 0; i < t->np; i++) {
        if (t->tag[i] != TAG_CONVEX) continue;
        if (band_size[i]) locked_size_sum += area[i];
        if (band_ang[i]) locked_ang_sum += area[i];
    }
    check("size_picks_bigger_faces", locked_size_sum > locked_ang_sum);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
