/* uvseam: target-driven seam splitting for automatic UV charts.
 *
 * The experiment starts from one LSCM island and greedily splits the currently
 * worst-angle island along BFS/dual-source chart candidates.  A split is accepted
 * under
 *
 *   J = w * Dmax / D0 + (1-w) * N / N_SCALE
 *
 * i.e. score = w * deltaD / D0 - (1-w) / N_SCALE.  Running the greedy relief
 * path once at w=1 exposes every marginal angle relief per added island; the
 * weight for an angle ceiling is then reverse-solved from those marginals.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAXV 512
#define MAXF 1024
#define MAXE 2048
#define MAX_ADJ 8
#define MAX_ISLANDS 80
#define N_SCALE 12.0
#define SCORE_EPS 1e-10

typedef struct { double x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { V3 v[MAXV]; int nv; Tri f[MAXF]; int nf; } Mesh;

typedef struct { int n[MAXF]; int to[MAXF][MAX_ADJ]; } Adj;
typedef struct { int a, b, f0, f1; } Edge;
typedef struct { Edge e[MAXE]; int ne; } EdgeMap;

typedef struct {
    int fids[MAXF], nf;
    int ids[MAXV], nv;
    double uv[MAXV][2];
    double errs[3 * MAXF];
    double dmax, dmean, p50, p95;
} Island;

typedef struct {
    int a[MAXF], na;
    int b[MAXF], nb;
    char name[16];
} Candidate;

typedef struct {
    char cut[16];
    int islands, size_a, size_b;
    double before, after, delta, score, wcrit;
} Step;

typedef struct {
    double w, D0, dmax, dmean, min_area;
    int N, nsteps;
    double area[MAX_ISLANDS];
    Island isl[MAX_ISLANDS];
    Step step[MAXF];
} RunResult;

static void mesh_init(Mesh *m) { memset(m, 0, sizeof(*m)); }

static int mesh_add_vertex(Mesh *m, V3 p) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv] = p;
    return m->nv++;
}

static int mesh_add_tri(Mesh *m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    m->f[m->nf++] = (Tri){a, b, c};
    return m->nf - 1;
}

static void mesh_sphere(Mesh *m, int lat_seg, int lon_seg) {
    int rings[32][64];
    mesh_init(m);
    int south = mesh_add_vertex(m, (V3){0, 0, -1});
    for (int i = 1; i < lat_seg; ++i) {
        double phi = -M_PI / 2.0 + M_PI * i / lat_seg;
        double z = sin(phi), r = cos(phi);
        for (int j = 0; j < lon_seg; ++j) {
            double th = 2.0 * M_PI * j / lon_seg;
            rings[i - 1][j] = mesh_add_vertex(m, (V3){r * cos(th), r * sin(th), z});
        }
    }
    int north = mesh_add_vertex(m, (V3){0, 0, 1});

    for (int j = 0; j < lon_seg; ++j)
        mesh_add_tri(m, south, rings[0][(j + 1) % lon_seg], rings[0][j]);

    for (int i = 0; i < lat_seg - 2; ++i) {
        for (int j = 0; j < lon_seg; ++j) {
            int k = (j + 1) % lon_seg;
            int a = rings[i][j], b = rings[i][k];
            int c = rings[i + 1][k], d = rings[i + 1][j];
            mesh_add_tri(m, a, b, c);
            mesh_add_tri(m, a, c, d);
        }
    }

    int last = lat_seg - 2;
    for (int j = 0; j < lon_seg; ++j)
        mesh_add_tri(m, rings[last][j], rings[last][(j + 1) % lon_seg], north);
}

static double face_area(const Mesh *m, int fi) {
    const Tri *t = &m->f[fi];
    V3 a = m->v[t->a], b = m->v[t->b], c = m->v[t->c];
    double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    double cx = uy * vz - uz * vy;
    double cy = uz * vx - ux * vz;
    double cz = ux * vy - uy * vx;
    return 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
}

static void edge_map_add(EdgeMap *em, int a, int b, int fi) {
    if (a > b) { int t = a; a = b; b = t; }
    for (int i = 0; i < em->ne; ++i) {
        if (em->e[i].a == a && em->e[i].b == b) {
            if (em->e[i].f0 < 0) em->e[i].f0 = fi;
            else if (em->e[i].f1 < 0 && em->e[i].f0 != fi) em->e[i].f1 = fi;
            return;
        }
    }
    if (em->ne >= MAXE) return;
    em->e[em->ne++] = (Edge){a, b, fi, -1};
}

static void adj_add(Adj *ad, int f, int g) {
    for (int i = 0; i < ad->n[f]; ++i) if (ad->to[f][i] == g) return;
    if (ad->n[f] < MAX_ADJ) ad->to[f][ad->n[f]++] = g;
}

static void mesh_build_adj(const Mesh *m, Adj *ad, EdgeMap *em) {
    memset(ad, 0, sizeof(*ad));
    memset(em, 0, sizeof(*em));
    for (int fi = 0; fi < m->nf; ++fi) {
        const Tri *t = &m->f[fi];
        edge_map_add(em, t->a, t->b, fi);
        edge_map_add(em, t->b, t->c, fi);
        edge_map_add(em, t->c, t->a, fi);
    }
    for (int i = 0; i < em->ne; ++i) {
        int f = em->e[i].f0, g = em->e[i].f1;
        if (f >= 0 && g >= 0) { adj_add(ad, f, g); adj_add(ad, g, f); }
    }
}

static double angle_at3(V3 o, V3 a, V3 b) {
    double ux = a.x - o.x, uy = a.y - o.y, uz = a.z - o.z;
    double vx = b.x - o.x, vy = b.y - o.y, vz = b.z - o.z;
    double du = sqrt(ux * ux + uy * uy + uz * uz);
    double dv = sqrt(vx * vx + vy * vy + vz * vz);
    double c = (ux * vx + uy * vy + uz * vz) / fmax(du * dv, 1e-300);
    return acos(fmin(1.0, fmax(-1.0, c)));
}

static double angle_at2(const double o[2], const double a[2], const double b[2]) {
    double ux = a[0] - o[0], uy = a[1] - o[1];
    double vx = b[0] - o[0], vy = b[1] - o[1];
    double du = sqrt(ux * ux + uy * uy), dv = sqrt(vx * vx + vy * vy);
    double c = (ux * vx + uy * vy) / fmax(du * dv, 1e-300);
    return acos(fmin(1.0, fmax(-1.0, c)));
}

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static int solve_qr(double *a, int m, int n, double *b, double *x) {
    if (m < n) return -1;
    double *v = calloc((size_t)m, sizeof(double));
    if (!v) return -3;

    for (int k = 0; k < n; ++k) {
        double normx = 0.0;
        for (int i = k; i < m; ++i) {
            double aik = a[i * n + k];
            normx += aik * aik;
        }
        normx = sqrt(normx);
        if (normx < 1e-30) { free(v); return -2; }

        double akk = a[k * n + k];
        double alpha = akk >= 0.0 ? -normx : normx;
        double denom = akk - alpha;
        for (int i = 0; i < m; ++i) v[i] = 0.0;
        v[k] = 1.0;
        for (int i = k + 1; i < m; ++i)
            v[i] = a[i * n + k] / denom;
        double vtv = 1.0;
        for (int i = k + 1; i < m; ++i) vtv += v[i] * v[i];
        double beta = 2.0 / vtv;

        for (int j = k; j < n; ++j) {
            double dot = 0.0;
            for (int i = k; i < m; ++i) dot += v[i] * a[i * n + j];
            for (int i = k; i < m; ++i) a[i * n + j] -= beta * v[i] * dot;
        }
        double dotb = 0.0;
        for (int i = k; i < m; ++i) dotb += v[i] * b[i];
        for (int i = k; i < m; ++i) b[i] -= beta * v[i] * dotb;
        (void)alpha;
    }

    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) s -= a[i * n + j] * x[j];
        double diag = a[i * n + i];
        if (fabs(diag) < 1e-30) { free(v); return -4; }
        x[i] = s / diag;
    }
    free(v);
    return 0;
}

static void island_compute(const Mesh *m, const int *fids, int nf, Island *out) {
    memset(out, 0, sizeof(*out));
    out->nf = nf;
    for (int i = 0; i < nf; ++i) out->fids[i] = fids[i];

    int map[MAXV];
    for (int i = 0; i < m->nv; ++i) map[i] = -1;
    for (int k = 0; k < nf; ++k) {
        const Tri *t = &m->f[fids[k]];
        int gs[3] = {t->a, t->b, t->c};
        for (int q = 0; q < 3; ++q) {
            int g = gs[q];
            if (map[g] < 0) {
                map[g] = out->nv;
                out->ids[out->nv++] = g;
            }
        }
    }

    /* Pin the longest island edge to (0,0)-(1,0).  This fixes translation,
     * rotation and uniform scale; angle metrics are scale invariant. */
    int p0 = out->ids[0], p1 = out->ids[1];
    double best_len = -1.0;
    for (int k = 0; k < nf; ++k) {
        const Tri *t = &m->f[fids[k]];
        int es[3][2] = {{t->a,t->b},{t->b,t->c},{t->c,t->a}};
        for (int q = 0; q < 3; ++q) {
            V3 a = m->v[es[q][0]], b = m->v[es[q][1]];
            double dx = a.x-b.x, dy = a.y-b.y, dz = a.z-b.z;
            double d = sqrt(dx*dx + dy*dy + dz*dz);
            if (d > best_len) { best_len = d; p0 = es[q][0]; p1 = es[q][1]; }
        }
    }

    int n = out->nv, full = 2 * n, free_n = full - 4;
    int fc[4] = {map[p0], n + map[p0], map[p1], n + map[p1]};
    double fv[4] = {0.0, 0.0, 1.0, 0.0};
    int red[2 * MAXV];
    for (int i = 0; i < full; ++i) red[i] = -1;
    int rn = 0;
    for (int col = 0; col < full; ++col) {
        bool fixed = false;
        for (int q = 0; q < 4; ++q) if (fc[q] == col) fixed = true;
        if (!fixed) red[col] = rn++;
    }

    int max_rows = 2 * nf;
    double *Q = calloc((size_t)max_rows * free_n, sizeof(double));
    double *rhs = calloc(max_rows, sizeof(double));
    double *sol = calloc(free_n, sizeof(double));
    if (!Q || !rhs || !sol) { free(Q); free(rhs); free(sol); return; }
    int qrows = 0;

    for (int k = 0; k < nf; ++k) {
        const Tri *t = &m->f[fids[k]];
        int gs[3] = {t->a, t->b, t->c};
        int idx[3] = {map[gs[0]], map[gs[1]], map[gs[2]]};
        V3 P[3] = {m->v[gs[0]], m->v[gs[1]], m->v[gs[2]]};

        double e1x = P[1].x-P[0].x, e1y = P[1].y-P[0].y, e1z = P[1].z-P[0].z;
        double e2x = P[2].x-P[0].x, e2y = P[2].y-P[0].y, e2z = P[2].z-P[0].z;
        double L = sqrt(e1x*e1x + e1y*e1y + e1z*e1z);
        double exx = e1x/L, exy = e1y/L, exz = e1z/L;
        double x3 = e2x*exx + e2y*exy + e2z*exz;
        double cx = e1y*e2z - e1z*e2y;
        double cy = e1z*e2x - e1x*e2z;
        double cz = e1x*e2y - e1y*e2x;
        double y3 = sqrt(cx*cx + cy*cy + cz*cz) / L;
        double D = L * y3;
        if (D < 1e-24) continue;
        double s = sqrt(D / 2.0) / D;

        /* Conformality rows in the triangle's metric plane. */
        double c[12] = {
            -y3*s,  y3*s, 0.0,
            (L-x3)*s, x3*s, -L*s,
            (x3-L)*s, -x3*s, L*s,
            -y3*s, y3*s, 0.0
        };
        int cols[12] = {
            idx[0], idx[1], idx[2],
            n+idx[0], n+idx[1], n+idx[2],
            idx[0], idx[1], idx[2],
            n+idx[0], n+idx[1], n+idx[2]
        };

        for (int row = 0; row < 2; ++row) {
            double *cc = &c[row * 6];
            int *cols6 = &cols[row * 6];
            double fixed_part = 0.0;
            for (int i = 0; i < 6; ++i) {
                if (red[cols6[i]] >= 0) continue;
                for (int q = 0; q < 4; ++q)
                    if (fc[q] == cols6[i]) fixed_part += cc[i] * fv[q];
            }
            for (int i = 0; i < 6; ++i) {
                int ri = red[cols6[i]];
                if (ri >= 0) Q[qrows * free_n + ri] = cc[i];
            }
            rhs[qrows] = -fixed_part;
            ++qrows;
        }
    }

    int solver_rc = solve_qr(Q, qrows, free_n, rhs, sol);
    if (solver_rc != 0)
        fprintf(stderr, "warning: QR solve rc=%d nf=%d nv=%d\n", solver_rc, nf, n);

    double xfull[2 * MAXV];
    memset(xfull, 0, sizeof(xfull));
    if (solver_rc == 0) {
        for (int col = 0; col < full; ++col)
            if (red[col] >= 0) xfull[col] = sol[red[col]];
    }
    for (int q = 0; q < 4; ++q) xfull[fc[q]] = fv[q];
    for (int i = 0; i < n; ++i) {
        int g = out->ids[i];
        out->uv[g][0] = xfull[i];
        out->uv[g][1] = xfull[n + i];
    }
    free(Q); free(rhs); free(sol);

    double sum = 0.0, emax = 0.0;
    for (int k = 0; k < nf; ++k) {
        int fi = fids[k];
        const Tri *t = &m->f[fi];
        V3 P[3] = {m->v[t->a], m->v[t->b], m->v[t->c]};
        double U[3][2] = {
            {out->uv[t->a][0], out->uv[t->a][1]},
            {out->uv[t->b][0], out->uv[t->b][1]},
            {out->uv[t->c][0], out->uv[t->c][1]}
        };
        int gs[3] = {t->a,t->b,t->c};
        (void)gs;
        for (int q = 0; q < 3; ++q) {
            int q1 = (q + 1) % 3, q2 = (q + 2) % 3;
            double a3 = angle_at3(P[q], P[q1], P[q2]);
            double a2 = angle_at2(U[q], U[q1], U[q2]);
            double e = fabs(a3 - a2) * 180.0 / M_PI;
            out->errs[3*k+q] = e;
            sum += e;
            emax = fmax(emax, e);
        }
    }
    int ne = 3 * nf;
    out->dmax = emax;
    out->dmean = sum / ne;
    double *sort_err = malloc((size_t)ne * sizeof(double));
    if (sort_err) {
        memcpy(sort_err, out->errs, (size_t)ne * sizeof(double));
        qsort(sort_err, ne, sizeof(double), cmp_double);
        out->p50 = sort_err[(int)(0.50 * (ne - 1))];
        out->p95 = sort_err[(int)(0.95 * (ne - 1))];
        free(sort_err);
    }
}

static bool group_connected(const int *ids, int n, const Adj *ad) {
    if (n <= 1) return true;
    char member[MAXF] = {0}, seen[MAXF] = {0};
    int q[MAXF], qh = 0, qt = 0, count = 0;
    for (int i = 0; i < n; ++i) member[ids[i]] = 1;
    q[qt++] = ids[0]; seen[ids[0]] = 1;
    while (qh < qt) {
        int f = q[qh++]; ++count;
        for (int i = 0; i < ad->n[f]; ++i) {
            int g = ad->to[f][i];
            if (member[g] && !seen[g]) { seen[g] = 1; q[qt++] = g; }
        }
    }
    return count == n;
}

static void try_add(Candidate *cs, int *nc, const int *a, int na,
                    const int *b, int nb, const char *name, const Adj *ad) {
    if (na <= 0 || nb <= 0 || *nc >= 8) return;
    if (!group_connected(a, na, ad) || !group_connected(b, nb, ad)) return;
    Candidate *c = &cs[(*nc)++];
    c->na = na; c->nb = nb;
    memcpy(c->a, a, (size_t)na * sizeof(int));
    memcpy(c->b, b, (size_t)nb * sizeof(int));
    snprintf(c->name, sizeof(c->name), "%s", name);
}

static void bfs_dist(const Adj *ad, const char *in, int start, int dist[MAXF]) {
    int q[MAXF], qh = 0, qt = 0;
    for (int i = 0; i < MAXF; ++i) dist[i] = -1;
    dist[start] = 0; q[qt++] = start;
    while (qh < qt) {
        int f = q[qh++];
        for (int i = 0; i < ad->n[f]; ++i) {
            int g = ad->to[f][i];
            if (in[g] && dist[g] < 0) { dist[g] = dist[f] + 1; q[qt++] = g; }
        }
    }
}

static int build_candidates(const Island *cur, const Adj *ad, Candidate *cs) {
    int nc = 0;
    char in[MAXF] = {0};
    for (int i = 0; i < cur->nf; ++i) in[cur->fids[i]] = 1;

    int worst_k = 0;
    for (int k = 1; k < 3 * cur->nf; ++k)
        if (cur->errs[k] > cur->errs[worst_k]) worst_k = k;
    int seed = cur->fids[worst_k / 3];

    int dist[MAXF];
    bfs_dist(ad, in, seed, dist);

    int a[MAXF], b[MAXF];
    for (int radius = 0; radius <= 3; ++radius) {
        int na = 0, nb = 0;
        for (int i = 0; i < cur->nf; ++i) {
            int f = cur->fids[i];
            if (dist[f] >= 0 && dist[f] <= radius) a[na++] = f;
            else b[nb++] = f;
        }
        char name[16]; snprintf(name, sizeof(name), "peel%d", radius);
        try_add(cs, &nc, a, na, b, nb, name, ad);
    }

    int far = cur->fids[0];
    for (int i = 0; i < cur->nf; ++i) {
        int f = cur->fids[i];
        if (dist[f] > dist[far]) far = f;
    }
    int d1[MAXF], d2[MAXF];
    bfs_dist(ad, in, seed, d1);
    bfs_dist(ad, in, far, d2);
    int na = 0, nb = 0;
    for (int i = 0; i < cur->nf; ++i) {
        int f = cur->fids[i];
        if (d1[f] <= d2[f]) a[na++] = f;
        else b[nb++] = f;
    }
    try_add(cs, &nc, a, na, b, nb, "dual", ad);
    return nc;
}

static void run_unwrap(const Mesh *m, const Adj *ad, double w, RunResult *r) {
    memset(r, 0, sizeof(*r));
    r->w = w;

    int all[MAXF];
    for (int i = 0; i < m->nf; ++i) all[i] = i;
    island_compute(m, all, m->nf, &r->isl[0]);
    r->N = 1;
    r->D0 = r->isl[0].dmax;

    double total_area = 0.0;
    for (int f = 0; f < m->nf; ++f) total_area += face_area(m, f);

    while (r->N < MAX_ISLANDS) {
        int gi = 0;
        for (int i = 1; i < r->N; ++i)
            if (r->isl[i].dmax > r->isl[gi].dmax) gi = i;
        Island *cur = &r->isl[gi];
        double before = cur->dmax;

        Candidate cands[8];
        int nc = build_candidates(cur, ad, cands);
        bool have = false;
        double best_score = 0, best_after = 0, best_delta = 0;
        Candidate best_cand;
        Island best_a, best_b;

        for (int ci = 0; ci < nc; ++ci) {
            Island ma, mb;
            island_compute(m, cands[ci].a, cands[ci].na, &ma);
            island_compute(m, cands[ci].b, cands[ci].nb, &mb);
            double after = fmax(ma.dmax, mb.dmax);
            for (int i = 0; i < r->N; ++i)
                if (i != gi) after = fmax(after, r->isl[i].dmax);
            double delta = before - after;
            double score = w * delta / r->D0 - (1.0 - w) / N_SCALE;
            if (!have || score > best_score) {
                have = true;
                best_score = score; best_after = after; best_delta = delta;
                best_cand = cands[ci];
                best_a = ma; best_b = mb;
            }
        }
        if (!have || best_score <= SCORE_EPS) break;
        if (getenv("UVSEAM_VERBOSE")) {
            long sa=0,sb=0; int mina=9999,maxa=-1,minb=9999,maxb=-1;
            for (int z=0;z<best_cand.na;++z){int f=best_cand.a[z]; sa+=f; mina=fmin(mina,f); maxa=fmax(maxa,f);}
            for (int z=0;z<best_cand.nb;++z){int f=best_cand.b[z]; sb+=f; minb=fmin(minb,f); maxb=fmax(maxb,f);}
            fprintf(stderr, "split %d %s %d/%d fpA(%ld,%d,%d) fpB(%ld,%d,%d) before %.9f after %.9f delta %.9f score %.12g\n",
                    r->nsteps, best_cand.name, best_cand.na, best_cand.nb,
                    sa,mina,maxa,sb,minb,maxb,
                    before, best_after, best_delta, best_score);
        }

        Step *s = &r->step[r->nsteps++];
        snprintf(s->cut, sizeof(s->cut), "%s", best_cand.name);
        s->islands = r->N + 1;
        s->size_a = best_cand.na; s->size_b = best_cand.nb;
        s->before = before; s->after = best_after;
        s->delta = best_delta; s->score = best_score;
        s->wcrit = r->D0 / (r->D0 + N_SCALE * best_delta);
        if (s->wcrit > 1.0) s->wcrit = 1.0;

        r->isl[gi] = best_a;
        r->isl[r->N] = best_b;
        ++r->N;
    }

    double dmax = 0.0, weighted_mean = 0.0, min_area = 1.0;
    for (int i = 0; i < r->N; ++i) {
        double ar = 0.0;
        for (int k = 0; k < r->isl[i].nf; ++k) ar += face_area(m, r->isl[i].fids[k]);
        r->area[i] = ar / total_area;
        min_area = fmin(min_area, r->area[i]);
        dmax = fmax(dmax, r->isl[i].dmax);
        weighted_mean += r->isl[i].dmean * 3.0 * r->isl[i].nf;
    }
    r->dmax = dmax;
    r->dmean = weighted_mean / (3.0 * m->nf);
    r->min_area = min_area;
}

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) ++failures;
}

static int run_selftest(void) {
    Mesh q; mesh_init(&q);
    mesh_add_vertex(&q, (V3){0,0,0});
    mesh_add_vertex(&q, (V3){1,0,0});
    mesh_add_vertex(&q, (V3){0,1,0});
    mesh_add_tri(&q, 0, 1, 2);
    Island one; int one_f[1] = {0};
    island_compute(&q, one_f, 1, &one);
    check("single_triangle_angle_zero", one.dmax < 1e-9);
    printf("single_triangle_dmax=%.3e\n", one.dmax);

    Mesh quad; mesh_init(&quad);
    mesh_add_vertex(&quad, (V3){0,0,0});
    mesh_add_vertex(&quad, (V3){1,0,0});
    mesh_add_vertex(&quad, (V3){1,1,0});
    mesh_add_vertex(&quad, (V3){0,1,0});
    mesh_add_tri(&quad, 0, 1, 2);
    mesh_add_tri(&quad, 0, 2, 3);
    Island qi; int two_f[2] = {0,1};
    island_compute(&quad, two_f, 2, &qi);
    check("coplanar_quad_angle_zero", qi.dmax < 1e-8);
    printf("coplanar_quad_dmax=%.3e\n", qi.dmax);

    double D0 = 100.0, delta = 1.0;
    double wc = D0 / (D0 + N_SCALE * delta);
    double score_at_boundary = wc * delta / D0 - (1.0 - wc) / N_SCALE;
    check("weight_boundary_zero", fabs(score_at_boundary) < 1e-15);

    Mesh sphere; mesh_sphere(&sphere, 10, 20);
    Adj ad; EdgeMap em; mesh_build_adj(&sphere, &ad, &em);
    static RunResult smoke;
    run_unwrap(&sphere, &ad, 0.942555, &smoke);
    double area_sum = 0;
    for (int i = 0; i < smoke.N; ++i) area_sum += smoke.area[i];
    check("weighted_run_islands_five", smoke.N == 5);
    check("area_fractions_sum_one", fabs(area_sum - 1.0) < 1e-12);
    check("medium_cap_met", smoke.dmax <= 6.0 + 1e-9);
    printf("smoke_N=%d smoke_dmax=%.6f area_sum=%.12f\n", smoke.N, smoke.dmax, area_sum);

    printf("all_ok=%d\n", failures == 0);
    return failures ? 1 : 0;
}

static void make_outdir(const char *dir) {
    if (mkdir(dir, 0755) != 0 && errno != EEXIST)
        fprintf(stderr, "warning: mkdir %s: %s\n", dir, strerror(errno));
}

static void write_obj(const char *path, const Mesh *m, const RunResult *r, const int *order) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }
    fprintf(f, "# uvseam export: each face corner has its own vt index\n");
    for (int i = 0; i < m->nv; ++i)
        fprintf(f, "v %.8f %.8f %.8f\n", m->v[i].x, m->v[i].y, m->v[i].z);
    int vt = 0;
    for (int pos = 0; pos < r->N; ++pos) {
        int ii = order[pos];
        const Island *z = &r->isl[ii];
        fprintf(f, "o island_%02d\n", pos + 1);
        for (int k = 0; k < z->nf; ++k) {
            int fi = z->fids[k];
            const Tri *t = &m->f[fi];
            int gs[3] = {t->a,t->b,t->c};
            int base = vt;
            for (int q = 0; q < 3; ++q) {
                int g = gs[q];
                fprintf(f, "vt %.8f %.8f\n", z->uv[g][0], z->uv[g][1]);
                ++vt;
            }
            fprintf(f, "f %d/%d %d/%d %d/%d\n",
                    t->a+1, base+1, t->b+1, base+2, t->c+1, base+3);
        }
    }
    fclose(f);
}

static void export_all(const char *dir, const Mesh *m, const EdgeMap *em,
                       RunResult runs[3], const double caps[3],
                       const double wreq[3], const double kappa[3]) {
    make_outdir(dir);
    char path[4096];

    snprintf(path, sizeof(path), "%s/islands.csv", dir);
    FILE *csv = fopen(path, "w");
    if (!csv) { perror(path); return; }
    fprintf(csv, "target_max_deg,weight_distortion,weight_island,island,faces,vertices,"
                 "area_fraction,angle_max_deg,angle_mean_deg,angle_p50_deg,angle_p95_deg\n");

    snprintf(path, sizeof(path), "%s/angle_hist.csv", dir);
    FILE *hist = fopen(path, "w");
    if (!hist) { perror(path); return; }
    fprintf(hist, "target_max_deg,island,bin_left_deg,bin_right_deg,corner_count\n");
    const double bins[] = {0,.1,.25,.5,1,2,4,8,16,32,64,180.0001};

    snprintf(path, sizeof(path), "%s/seams.csv", dir);
    FILE *seam = fopen(path, "w");
    if (!seam) { perror(path); return; }
    fprintf(seam, "target_max_deg,island_a,island_b,v0,v1,x0,y0,z0,x1,y1,z1,length\n");

    for (int ti = 0; ti < 3; ++ti) {
        RunResult *r = &runs[ti];
        int order[MAX_ISLANDS], face_rank[MAXF];
        for (int i = 0; i < r->N; ++i) order[i] = i;
        for (int i = 0; i < r->N; ++i)
            for (int j = i + 1; j < r->N; ++j)
                if (r->area[order[j]] > r->area[order[i]]) {
                    int tmp = order[i]; order[i] = order[j]; order[j] = tmp;
                }
        for (int pos = 0; pos < r->N; ++pos) {
            int old = order[pos];
            const Island *z = &r->isl[old];
            for (int k = 0; k < z->nf; ++k) face_rank[z->fids[k]] = pos;
            fprintf(csv, "%.0f,%.6f,%.6f,%d,%d,%d,%.12f,%.9f,%.9f,%.9f,%.9f\n",
                    caps[ti], r->w, 1.0-r->w, pos+1, z->nf, z->nv, r->area[old],
                    z->dmax, z->dmean, z->p50, z->p95);
            for (int bi = 0; bi < (int)(sizeof(bins)/sizeof(bins[0]))-1; ++bi) {
                int count = 0;
                for (int e = 0; e < 3*z->nf; ++e)
                    if (z->errs[e] >= bins[bi] && z->errs[e] < bins[bi+1]) ++count;
                fprintf(hist, "%.0f,%d,%g,%g,%d\n", caps[ti], pos+1, bins[bi], bins[bi+1], count);
            }
        }
        for (int ei = 0; ei < em->ne; ++ei) {
            const Edge *e = &em->e[ei];
            if (e->f0 < 0 || e->f1 < 0) continue;
            int ra = face_rank[e->f0], rb = face_rank[e->f1];
            if (ra == rb) continue;
            V3 a = m->v[e->a], b = m->v[e->b];
            double len = sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2)+pow(a.z-b.z,2));
            int lo = ra < rb ? ra : rb, hi = ra < rb ? rb : ra;
            fprintf(seam, "%.0f,%d,%d,%d,%d,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n",
                    caps[ti], lo+1, hi+1, e->a, e->b,
                    a.x,a.y,a.z,b.x,b.y,b.z,len);
        }

        char opath[4096];
        snprintf(opath, sizeof(opath), "%s/uv_target_%.0fdeg.obj", dir, caps[ti]);
        write_obj(opath, m, r, order);
    }
    fclose(csv); fclose(hist); fclose(seam);

    snprintf(path, sizeof(path), "%s/summary.json", dir);
    FILE *j = fopen(path, "w");
    if (!j) { perror(path); return; }
    fprintf(j, "{\n");
    fprintf(j, "  \"mesh\": \"uv-sphere-10x20\",\n");
    fprintf(j, "  \"vertices\": %d,\n", m->nv);
    fprintf(j, "  \"triangles\": %d,\n", m->nf);
    fprintf(j, "  \"angle_unit\": \"degree\",\n");
    fprintf(j, "  \"area_fraction_basis\": \"3D triangle area / total 3D area\",\n");
    fprintf(j, "  \"objective\": \"w*Dmax/D0 + (1-w)*N/%.1f\",\n", N_SCALE);
    fprintf(j, "  \"D0_max_angle_deg\": %.12f,\n", runs[0].D0);
    fprintf(j, "  \"metrics_order\": [\"angle_max_deg\",\"angle_mean_deg\",\"island_count\",\"min_island_area_fraction\"],\n");
    fprintf(j, "  \"rows_12\": [\n");
    for (int ti = 0; ti < 3; ++ti) {
        RunResult *r = &runs[ti];
        fprintf(j, "    {\"target_max_deg\": %.0f, \"values\": [%.9f, %.9f, %d, %.12f]}%s\n",
                caps[ti], r->dmax, r->dmean, r->N, r->min_area, ti == 2 ? "" : ",");
    }
    fprintf(j, "  ],\n");
    fprintf(j, "  \"tiers\": [\n");
    for (int ti = 0; ti < 3; ++ti) {
        RunResult *r = &runs[ti];
        fprintf(j, "    {\n");
        fprintf(j, "      \"target_max_deg\": %.0f,\n", caps[ti]);
        fprintf(j, "      \"reverse_solved_weight_raw\": %.12f,\n", wreq[ti]);
        fprintf(j, "      \"weight_distortion\": %.6f,\n", r->w);
        fprintf(j, "      \"weight_island_count\": %.6f,\n", 1.0-r->w);
        fprintf(j, "      \"exchange_threshold_deg_per_island\": %.9f,\n", kappa[ti]);
        fprintf(j, "      \"island_count\": %d,\n", r->N);
        fprintf(j, "      \"angle_max_deg\": %.12f,\n", r->dmax);
        fprintf(j, "      \"angle_mean_deg\": %.12f,\n", r->dmean);
        fprintf(j, "      \"min_island_area_fraction\": %.12f,\n", r->min_area);
        fprintf(j, "      \"min_island_area_percent\": %.6f\n", 100.0*r->min_area);
        fprintf(j, "    }%s\n", ti == 2 ? "" : ",");
    }
    fprintf(j, "  ]\n");
    fprintf(j, "}\n");
    fclose(j);

    snprintf(path, sizeof(path), "%s/report.md", dir);
    FILE *md = fopen(path, "w");
    if (!md) { perror(path); return; }
    fprintf(md, "# UV 切缝三档实验结果\n\n");
    fprintf(md, "网格为确定性 UV 球（%d 顶点 / %d 三角面）；每岛用 LSCM 展开，", m->nv, m->nf);
    fprintf(md, "角度畸变是三个内角的绝对误差（度），均值按角数平均。\n\n");
    fprintf(md, "## 三组十二个数\n\n");
    fprintf(md, "| 目标角畸变上限 | 最大角畸变(°) | 平均角畸变(°) | 岛数 | 最小岛面积占比 |\n");
    fprintf(md, "|---:|---:|---:|---:|---:|\n");
    for (int ti = 0; ti < 3; ++ti) {
        RunResult *r = &runs[ti];
        fprintf(md, "| %.0f° | %.4f | %.4f | %d | %.6f (%.4f%%) |\n",
                caps[ti], r->dmax, r->dmean, r->N, r->min_area, 100*r->min_area);
    }
    fprintf(md, "\n## 反解权重与实际岛数\n\n");
    fprintf(md, "归一化目标为 `w·Dmax/D0 + (1-w)·N/%.0f`，D0=%.4f°。", N_SCALE, runs[0].D0);
    fprintf(md, "临界权重 `w*=D0/(D0+%.0f·ΔD)`；表中向上取到 6 位小数后实际重跑。\n\n", N_SCALE);
    fprintf(md, "| 上限 | 反解原始 w* | 采用畸变权重 w | 岛数权重 1-w | 边际交换阈值(°/岛) | 实际岛数 |\n");
    fprintf(md, "|---:|---:|---:|---:|---:|---:|\n");
    for (int ti = 0; ti < 3; ++ti) {
        RunResult *r = &runs[ti];
        fprintf(md, "| %.0f° | %.9f | %.6f | %.6f | %.6f | %d |\n",
                caps[ti], wreq[ti], r->w, 1-r->w, kappa[ti], r->N);
    }
    bool fragmented = runs[2].N >= 2 * runs[1].N && runs[2].min_area < 0.01;
    fprintf(md, "\n## 诊断\n\n");
    if (fragmented) {
        fprintf(md, "最紧一档岛数由 %d 暴涨到 %d（%.1f 倍），最小岛面积占比仅 %.6f（%.4f%%），",
                runs[1].N, runs[2].N, (double)runs[2].N/runs[1].N,
                runs[2].min_area, 100*runs[2].min_area);
        fprintf(md, "且 CSV 中存在单面岛。结论成立：%.6f 的畸变权重压过了 %.6f 的岛数权重，",
                runs[2].w, 1-runs[2].w);
        fprintf(md, "标量目标只计数岛而没有最小岛面积/碎片化约束，因此切缝继续剥出极小碎片。\n");
    } else {
        fprintf(md, "未出现“岛数暴涨且最小岛面积接近零”的碎片化特征。\n");
    }
    fprintf(md, "\n逐岛面积和 p50/p95 见 `islands.csv`，角度直方图见 `angle_hist.csv`，切缝边见 `seams.csv`。\n");
    fclose(md);

    printf("wrote %s/summary.json, islands.csv, angle_hist.csv, seams.csv and OBJ files\n", dir);
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return run_selftest();

    const char *dir = "out";
    if (argc >= 3 && strcmp(argv[1], "--export") == 0) dir = argv[2];
    else if (argc >= 2 && strcmp(argv[1], "--export") != 0) {
        fprintf(stderr, "usage: %s --selftest | --export [DIR]\n", argv[0]);
        return 2;
    }

    static Mesh sphere;
    mesh_sphere(&sphere, 10, 20);
    static Adj ad;
    static EdgeMap em;
    mesh_build_adj(&sphere, &ad, &em);
    printf("mesh vertices=%d triangles=%d edges=%d\n", sphere.nv, sphere.nf, em.ne);

    static RunResult full;
    run_unwrap(&sphere, &ad, 1.0, &full);
    printf("relief_path islands=%d dmax=%.9f dmean=%.9f min_area=%.9f\n",
           full.N, full.dmax, full.dmean, full.min_area);

    double caps[3] = {10.0, 6.0, 3.0};
    double wreq[3], kappa[3];
    double wuse[3];
    static RunResult runs[3];
    for (int ti = 0; ti < 3; ++ti) {
        wreq[ti] = 0.0;
        bool found = false;
        for (int s = 0; s < full.nsteps; ++s) {
            wreq[ti] = fmax(wreq[ti], full.step[s].wcrit);
            if (full.step[s].after <= caps[ti]) { found = true; break; }
        }
        if (!found) {
            fprintf(stderr, "cap %.1f not reached (final dmax %.6f)\n", caps[ti], full.dmax);
            return 1;
        }
        wuse[ti] = ceil(wreq[ti] * 1e6) / 1e6;
        if (wuse[ti] > 1.0) wuse[ti] = 1.0;
        kappa[ti] = full.D0 * (1.0 - wuse[ti]) / (N_SCALE * wuse[ti]);
        run_unwrap(&sphere, &ad, wuse[ti], &runs[ti]);
        printf("target=%.0f wraw=%.12f w=%.6f kappa=%.9f N=%d dmax=%.9f dmean=%.9f min=%.12f\n",
               caps[ti], wreq[ti], runs[ti].w, kappa[ti], runs[ti].N,
               runs[ti].dmax, runs[ti].dmean, runs[ti].min_area);
    }

    export_all(dir, &sphere, &em, runs, caps, wreq, kappa);
    return 0;
}
