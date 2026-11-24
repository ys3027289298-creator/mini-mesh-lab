
/* triangulate: concave polygon ear clipping with an aspect-ratio gate.
 *
 * The gate is derived from the polygon's mean edge length L and a target
 * minimum triangle area A*.  For a triangle whose longest edge is l, the
 * altitude to that edge is h = 2A/l, hence
 *
 *     aspect = l/h = l^2 / (2A).
 *
 * Using the polygon-wide mean L in place of the local l gives the requested
 * global threshold:
 *
 *     gate = L^2 / (2A*).
 *
 * If no convex, containment-free ear passes that gate, the current loop is
 * split along a balanced internal diagonal instead of clipping a sliver.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#define MAXV 64
#define TARGET_MIN_AREA 0.01
#define DEGEN_AREA 1.0e-8
#define GEOM_EPS 1.0e-12
#define ON_EPS 1.0e-10

typedef struct { double x, y; } Vec2;

typedef struct {
    int idx[3];
    double edge[3];
    double longest;
    double area;
    double aspect;
    int gated_ear;   /* 1 = accepted by the ear gate, 0 = terminal triangle */
} Tri;

typedef struct {
    Vec2 p[MAXV];
    int n;
    Tri tris[MAXV - 2];
    int ntri;
    int splits;
    double mean_edge;
    double gate;
    double target_area;
    double degen_area;
    double min_area;
    int degenerate_count;
    int first_degenerate;
} TriReport;

static double cross3(Vec2 a, Vec2 b, Vec2 c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static double edge_len(Vec2 a, Vec2 b) {
    return hypot(a.x - b.x, a.y - b.y);
}

static double signed_area(const Vec2* p, int n) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        Vec2 a = p[i], b = p[(i + 1) % n];
        s += a.x * b.y - a.y * b.x;
    }
    return 0.5 * s;
}

static int on_segment(Vec2 q, Vec2 a, Vec2 b) {
    if (fabs(cross3(a, b, q)) > ON_EPS) return 0;
    return q.x >= fmin(a.x, b.x) - ON_EPS &&
           q.x <= fmax(a.x, b.x) + ON_EPS &&
           q.y >= fmin(a.y, b.y) - ON_EPS &&
           q.y <= fmax(a.y, b.y) + ON_EPS;
}

/* A point is inside or on the triangle.  Containment includes the boundary:
 * another loop vertex on an ear diagonal makes that ear unsafe to clip.
 */
static int point_in_triangle(Vec2 q, Vec2 a, Vec2 b, Vec2 c) {
    double c1 = cross3(a, b, q);
    double c2 = cross3(b, c, q);
    double c3 = cross3(c, a, q);
    int neg = (c1 < -GEOM_EPS) || (c2 < -GEOM_EPS) || (c3 < -GEOM_EPS);
    int pos = (c1 > GEOM_EPS) || (c2 > GEOM_EPS) || (c3 > GEOM_EPS);
    return !(neg && pos);
}

static int point_on_loop(Vec2 q, const int* loop, int m, const Vec2* p) {
    for (int i = 0; i < m; ++i) {
        if (on_segment(q, p[loop[i]], p[loop[(i + 1) % m]])) return 1;
    }
    return 0;
}

static int point_in_loop(Vec2 q, const int* loop, int m, const Vec2* p) {
    int inside = 0;
    for (int i = 0; i < m; ++i) {
        Vec2 a = p[loop[i]];
        Vec2 b = p[loop[(i + 1) % m]];
        if ((a.y > q.y) != (b.y > q.y)) {
            double xhit = a.x + (b.x - a.x) * (q.y - a.y) / (b.y - a.y);
            if (xhit > q.x) inside = !inside;
        }
    }
    return inside;
}

static int segments_intersect(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
    double c1 = cross3(a, b, c);
    double c2 = cross3(a, b, d);
    double c3 = cross3(c, d, a);
    double c4 = cross3(c, d, b);

    if (((c1 > GEOM_EPS && c2 < -GEOM_EPS) ||
         (c1 < -GEOM_EPS && c2 > GEOM_EPS)) &&
        ((c3 > GEOM_EPS && c4 < -GEOM_EPS) ||
         (c3 < -GEOM_EPS && c4 > GEOM_EPS))) {
        return 1;
    }
    return on_segment(c, a, b) || on_segment(d, a, b) ||
           on_segment(a, c, d) || on_segment(b, c, d);
}

static int diagonal_valid(const TriReport* r, const int* loop, int m,
                          int ipos, int jpos) {
    int ia = loop[ipos], ib = loop[jpos];
    Vec2 a = r->p[ia], b = r->p[ib];
    if (edge_len(a, b) <= GEOM_EPS) return 0;

    for (int k = 0; k < m; ++k) {
        int u = loop[k];
        int v = loop[(k + 1) % m];
        /* Edges incident to either diagonal endpoint meet it legally. */
        if (u == ia || v == ia || u == ib || v == ib) continue;
        if (segments_intersect(a, b, r->p[u], r->p[v])) return 0;
    }

    Vec2 mid = { 0.5 * (a.x + b.x), 0.5 * (a.y + b.y) };
    if (point_on_loop(mid, loop, m, r->p)) return 0;
    return point_in_loop(mid, loop, m, r->p);
}

static void triangle_quality(Tri* t, int ia, int ib, int ic,
                             const Vec2* p, int gated_ear) {
    t->idx[0] = ia;
    t->idx[1] = ib;
    t->idx[2] = ic;
    t->edge[0] = edge_len(p[ia], p[ib]);
    t->edge[1] = edge_len(p[ib], p[ic]);
    t->edge[2] = edge_len(p[ic], p[ia]);
    t->longest = fmax(t->edge[0], fmax(t->edge[1], t->edge[2]));
    t->area = fabs(cross3(p[ia], p[ib], p[ic])) * 0.5;
    t->aspect = (t->area > 1.0e-14)
        ? t->longest * t->longest / (2.0 * t->area)
        : INFINITY;
    t->gated_ear = gated_ear;
}

static void report_add_triangle(TriReport* r, int ia, int ib, int ic,
                                int gated_ear) {
    Tri* t = &r->tris[r->ntri];
    triangle_quality(t, ia, ib, ic, r->p, gated_ear);
    if (r->ntri == 0 || t->area < r->min_area) r->min_area = t->area;
    if (t->area < r->degen_area) {
        r->degenerate_count++;
        if (r->first_degenerate < 0) r->first_degenerate = r->ntri;
    }
    r->ntri++;
}

static int triangulate_rec(TriReport* r, int* loop, int m) {
    if (m == 3) {
        report_add_triangle(r, loop[0], loop[1], loop[2], 0);
        return 0;
    }
    if (m < 3 || r->ntri >= MAXV - 2) return -1;

    /* Pick the lowest-aspect convex ear that contains no remaining vertex and
     * passes the polygon-wide aspect gate.
     */
    int best_q = -1;
    double best_aspect = INFINITY;

    for (int q = 0; q < m; ++q) {
        int ip = loop[(q - 1 + m) % m];
        int iq = loop[q];
        int in = loop[(q + 1) % m];
        Vec2 a = r->p[ip], b = r->p[iq], c = r->p[in];

        if (cross3(a, b, c) <= GEOM_EPS) continue;

        int blocked = 0;
        for (int k = 0; k < m; ++k) {
            int id = loop[k];
            if (id == ip || id == iq || id == in) continue;
            if (point_in_triangle(r->p[id], a, b, c)) {
                blocked = 1;
                break;
            }
        }
        if (blocked) continue;

        Tri candidate;
        triangle_quality(&candidate, ip, iq, in, r->p, 1);
        if (candidate.aspect <= r->gate && candidate.aspect < best_aspect) {
            best_aspect = candidate.aspect;
            best_q = q;
        }
    }

    if (best_q >= 0) {
        int ip = loop[(best_q - 1 + m) % m];
        int iq = loop[best_q];
        int in = loop[(best_q + 1) % m];
        report_add_triangle(r, ip, iq, in, 1);

        int next[MAXV];
        int nm = 0;
        for (int q = 0; q < m; ++q) {
            if (q != best_q) next[nm++] = loop[q];
        }
        return triangulate_rec(r, next, nm);
    }

    /* No acceptable ear: split along the most balanced internal diagonal. */
    int bi = -1, bj = -1;
    int best_balance = INT_MAX;
    double best_cap_area = -1.0;
    for (int i = 0; i < m; ++i) {
        for (int j = i + 2; j < m; ++j) {
            if (i == 0 && j == m - 1) continue; /* already a boundary edge */
            if (!diagonal_valid(r, loop, m, i, j)) continue;

            int left_n = j - i + 1;
            int right_n = m - left_n + 2;
            int balance = abs(left_n - right_n);

            /* Among equally balanced cuts, avoid the diagonal that merely
             * isolates a microscopic cap triangle.
             */
            int im = (i - 1 + m) % m;
            int ip = (i + 1) % m;
            int jm = (j - 1 + m) % m;
            int jp = (j + 1) % m;
            double cap_area = fmin(
                fabs(cross3(r->p[loop[i]], r->p[loop[ip]], r->p[loop[j]])),
                fabs(cross3(r->p[loop[jm]], r->p[loop[j]], r->p[loop[i]])));
            cap_area = fmin(cap_area,
                fabs(cross3(r->p[loop[j]], r->p[loop[jp]], r->p[loop[i]])));
            cap_area = fmin(cap_area,
                fabs(cross3(r->p[loop[im]], r->p[loop[i]], r->p[loop[j]]))) * 0.5;

            if (balance < best_balance ||
                (balance == best_balance && cap_area > best_cap_area)) {
                best_balance = balance;
                best_cap_area = cap_area;
                bi = i;
                bj = j;
            }
        }
    }

    if (bi < 0) return -1;
    r->splits++;

    int left[MAXV], right[MAXV];
    int ln = 0, rn = 0;
    for (int k = bi; k <= bj; ++k) left[ln++] = loop[k];
    for (int k = bj; k < m; ++k) right[rn++] = loop[k];
    for (int k = 0; k <= bi; ++k) right[rn++] = loop[k];

    if (triangulate_rec(r, left, ln) != 0) return -1;
    return triangulate_rec(r, right, rn);
}

static int triangulate_polygon(const Vec2* p, int n,
                               double target_area, double degen_area,
                               TriReport* r) {
    if (n < 3 || n > MAXV || target_area <= 0.0 || degen_area <= 0.0) return -1;
    memset(r, 0, sizeof(*r));
    r->n = n;
    r->target_area = target_area;
    r->degen_area = degen_area;
    r->first_degenerate = -1;
    r->min_area = INFINITY;
    memcpy(r->p, p, sizeof(Vec2) * (size_t)n);

    if (signed_area(r->p, n) < 0.0) {
        for (int i = 0, j = n - 1; i < j; ++i, --j) {
            Vec2 tmp = r->p[i];
            r->p[i] = r->p[j];
            r->p[j] = tmp;
        }
    }

    double perimeter = 0.0;
    for (int i = 0; i < n; ++i) {
        perimeter += edge_len(r->p[i], r->p[(i + 1) % n]);
    }
    r->mean_edge = perimeter / (double)n;
    r->gate = r->mean_edge * r->mean_edge / (2.0 * target_area);

    int loop[MAXV];
    for (int i = 0; i < n; ++i) loop[i] = i;
    return triangulate_rec(r, loop, n);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static Vec2 compact_polygon[7];
static Vec2 medium_polygon[7];
static Vec2 slender_polygon[7];

static void make_case(Vec2* out, double width) {
    const double h = 1.0;
    const double ell = 2.0e-4;
    const double kink_aspect = 1000.0;
    const double bump_h = ell / kink_aspect;
    double x0 = 0.5 * (width - ell);

    out[0] = (Vec2){ 0.0, 0.0 };
    out[1] = (Vec2){ width, 0.0 };
    out[2] = (Vec2){ width, h };
    out[3] = (Vec2){ x0 + ell, h };
    out[4] = (Vec2){ x0 + 0.5 * ell, h + bump_h };
    out[5] = (Vec2){ x0, h };
    out[6] = (Vec2){ 0.0, h };
}

static void print_case(const char* name, const TriReport* r) {
    printf("\n[%s]\n", name);
    printf("%s_vertices=%d\n", name, r->n);
    printf("%s_mean_edge=%.12g\n", name, r->mean_edge);
    printf("%s_gate=%.12g\n", name, r->gate);
    printf("%s_triangle_count=%d\n", name, r->ntri);
    printf("%s_min_area=%.12g\n", name, r->min_area);
    printf("%s_degenerate_count=%d\n", name, r->degenerate_count);
    printf("%s_split_count=%d\n", name, r->splits);
}

int main(void) {
    make_case(compact_polygon, 2.0);
    make_case(medium_polygon, 8.0);
    make_case(slender_polygon, 32.0);

    TriReport rc, rm, rs;
    int okc = triangulate_polygon(compact_polygon, 7,
                                  TARGET_MIN_AREA, DEGEN_AREA, &rc);
    int okm = triangulate_polygon(medium_polygon, 7,
                                  TARGET_MIN_AREA, DEGEN_AREA, &rm);
    int oks = triangulate_polygon(slender_polygon, 7,
                                  TARGET_MIN_AREA, DEGEN_AREA, &rs);

    check("compact_triangulated", okc == 0);
    check("medium_triangulated", okm == 0);
    check("slender_triangulated", oks == 0);

    print_case("compact", &rc);
    print_case("medium", &rm);
    print_case("slender", &rs);

    printf("\nnine_numbers={%d,%.12g,%d; %d,%.12g,%d; %d,%.12g,%d}\n",
           rc.ntri, rc.min_area, rc.degenerate_count,
           rm.ntri, rm.min_area, rm.degenerate_count,
           rs.ntri, rs.min_area, rs.degenerate_count);

    /* Conditional diagnostic: global mean edge versus the local ear edges. */
    if (rs.degenerate_count > 0 && rs.first_degenerate >= 0) {
        const Tri* bad = &rs.tris[rs.first_degenerate];
        double local_gate = bad->longest * bad->longest /
                            (2.0 * TARGET_MIN_AREA);
        printf("\nslender_local_vs_gate:\n");
        printf("local_edges=%.12g,%.12g,%.12g\n",
               bad->edge[0], bad->edge[1], bad->edge[2]);
        printf("local_longest=%.12g\n", bad->longest);
        printf("bad_aspect=%.12g\n", bad->aspect);
        printf("global_mean_edge=%.12g\n", rs.mean_edge);
        printf("global_gate=%.12g\n", rs.gate);
        printf("local_gate_from_longest=%.12g\n", local_gate);
        printf("bad_aspect_passes_global_gate=%d\n",
               bad->aspect <= rs.gate ? 1 : 0);
    }

    check("triangle_count_compact", rc.ntri == 5);
    check("triangle_count_medium", rm.ntri == 5);
    check("triangle_count_slender", rs.ntri == 5);
    check("compact_zero_degenerate", rc.degenerate_count == 0);
    check("medium_zero_degenerate", rm.degenerate_count == 0);
    check("slender_global_gate_reproduces_one_degenerate",
          rs.degenerate_count == 1);
    check("compact_used_split_fallback", rc.splits >= 1);
    check("medium_used_split_fallback", rm.splits >= 1);

    if (oks == 0 && rs.first_degenerate >= 0) {
        const Tri* bad = &rs.tris[rs.first_degenerate];
        check("slender_bad_ear_passes_global_gate", bad->aspect <= rs.gate);
        check("global_gate_much_looser_than_local",
              rs.gate > 1000.0 * bad->longest * bad->longest /
                            (2.0 * TARGET_MIN_AREA));
    } else {
        check("slender_bad_ear_passes_global_gate", 0);
        check("global_gate_much_looser_than_local", 0);
    }

    /* Every triangulation of an n-vertex simple polygon has n-2 triangles and
     * preserves polygon area.
     */
    const struct { const char* name; const Vec2* p; const TriReport* r; } cases[] = {
        { "compact", compact_polygon, &rc },
        { "medium", medium_polygon, &rm },
        { "slender", slender_polygon, &rs },
    };
    for (int k = 0; k < 3; ++k) {
        const TriReport* r = cases[k].r;
        double sum = 0.0;
        int index_ok = 1;
        for (int i = 0; i < r->ntri; ++i) {
            sum += r->tris[i].area;
            for (int a = 0; a < 3; ++a) {
                int id = r->tris[i].idx[a];
                if (id < 0 || id >= r->n) index_ok = 0;
                for (int b = a + 1; b < 3; ++b) {
                    if (id == r->tris[i].idx[b]) index_ok = 0;
                }
            }
        }
        double pa = fabs(signed_area(cases[k].p, r->n));
        double tol = 1.0e-8 * fmax(1.0, pa);
        char label[64];
        snprintf(label, sizeof(label), "%s_indices_valid", cases[k].name);
        check(label, index_ok);
        snprintf(label, sizeof(label), "%s_area_preserved", cases[k].name);
        check(label, fabs(sum - pa) <= tol);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
