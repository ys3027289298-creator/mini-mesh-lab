/* falloff: soft-selection brush falloff vs mesh density.
 *
 * Brush of radius R at the origin, hand-entered exponent p:
 *
 *     w(r) = 1 - (r/R)^p   for r < R,   0 outside
 *     w(0) = 1: the centre is fully weighted, as intended.
 *
 * Both conditions on p are read off a square lattice of spacing h.  Let k = h/R.
 *
 *   A  vertex count.  A vertex moves when w >= W_MOVE, i.e. inside
 *      r_cut = R * W_MOVE^(1/p).  Needing N_MIN vertices asks r_cut to reach the
 *      first axial ring at r = h:
 *          W_MOVE^(1/p) >= k    =>    p >= p_A = ln(W_MOVE) / ln(k)
 *
 *   B  no visible step.  Adjacent weights may differ by at most G_MAX.  Two edges
 *      carry the largest jumps, one at each end of the falloff:
 *          centre edge  w(0)-w(h)   = k^p          <= G_MAX  =>  p >= p_CORE
 *          rim edge     w(R-h)-w(R) = 1-(1-k)^p    <= G_MAX  =>  p <= p_RIM
 *      with
 *          p_CORE = ln(G_MAX) / ln(k),
 *          p_RIM  = ln(1-G_MAX) / ln(1-k).
 *      B is the bracket [p_CORE, p_RIM], non-empty iff k < G_MAX.
 *
 * Usable exponent = [max(p_A, p_CORE), p_RIM].
 *
 * R=1, W_MOVE=0.5, N_MIN=5 (centre + first ring), G_MAX=0.3.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define R       1.0
#define W_MOVE  0.5
#define N_MIN   5
#define G_MAX   0.3
#define P_FIXED 2.0      /* the hand-entered exponent */
#define EXT     3.0

typedef struct { double h; int n; double* w; } Lattice;

static int lat_alloc(Lattice* L, double h) {
    L->h = h;
    L->n = (int)llround(2.0 * EXT / h) + 1;
    L->w = malloc(sizeof(double) * (size_t)L->n * (size_t)L->n);
    return L->w ? 0 : -1;
}

static void lat_free(Lattice* L) { free(L->w); L->w = NULL; }

static void lat_fill(Lattice* L, double p) {
    int n = L->n;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double x = -EXT + i * L->h;
            double y = -EXT + j * L->h;
            double r = sqrt(x*x + y*y);
            L->w[j*n + i] = (r >= R) ? 0.0 : 1.0 - pow(r / R, p);
        }
    }
}

static int affected(const Lattice* L) {
    int n = L->n, c = 0;
    for (int q = 0; q < n*n; q++) if (L->w[q] >= W_MOVE - 1e-12) c++;
    return c;
}

/* max adjacent-vertex jump and the midpoint radius of that edge */
static double edge_gradient(const Lattice* L, double* r_at) {
    int n = L->n;
    double g = 0.0, rg = 0.0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double w = L->w[j*n + i];
            double x = -EXT + i * L->h, y = -EXT + j * L->h;
            if (i + 1 < n) {
                double d = fabs(w - L->w[j*n + i + 1]);
                double r2 = sqrt((x + L->h)*(x + L->h) + y*y);
                if (d > g) { g = d; rg = 0.5 * (sqrt(x*x + y*y) + r2); }
            }
            if (j + 1 < n) {
                double d = fabs(w - L->w[(j+1)*n + i]);
                double r2 = sqrt(x*x + (y + L->h)*(y + L->h));
                if (d > g) { g = d; rg = 0.5 * (sqrt(x*x + y*y) + r2); }
            }
        }
    }
    if (r_at) *r_at = rg;
    return g;
}

/* cells whose corner jump exceeds G_MAX: the visible stair zone */
static double seam_area(const Lattice* L, int* cells_out) {
    int n = L->n, cells = 0;
    for (int j = 0; j + 1 < n; j++) {
        for (int i = 0; i + 1 < n; i++) {
            double a = L->w[j*n + i],     b = L->w[j*n + i + 1];
            double c = L->w[(j+1)*n + i], d = L->w[(j+1)*n + i + 1];
            double m = fabs(a-b);
            if (fabs(c-d) > m) m = fabs(c-d);
            if (fabs(a-c) > m) m = fabs(a-c);
            if (fabs(b-d) > m) m = fabs(b-d);
            if (m > G_MAX) cells++;
        }
    }
    if (cells_out) *cells_out = cells;
    return cells * L->h * L->h;
}

/* smallest p with N(p) >= N_MIN (N rises with p) */
static double p_coverage(Lattice* L) {
    double lo = 1e-4, hi = 64.0;
    for (int it = 0; it < 100; it++) {
        double mid = 0.5 * (lo + hi);
        lat_fill(L, mid);
        if (affected(L) >= N_MIN) hi = mid; else lo = mid;
    }
    return hi;
}

static double dmax(double a, double b) { return a > b ? a : b; }

/* the two edges that carry the largest jumps */
static double jump_core(double k, double p) { return pow(k, p); }
static double jump_rim(double k, double p)  { return 1.0 - pow(1.0 - k, p); }

int main(void) {
    static const double hs[] = { 0.5, 0.25, 0.125 };
    static const char* names[] = { "coarse", "medium", "fine" };
    enum { ND = 3 };
    int failures = 0;

    printf("falloff: R=%.3f W_MOVE=%.3f N_MIN=%d G_MAX=%.3f p=%.3f\n", R, W_MOVE, N_MIN, G_MAX, P_FIXED);
    printf("k=h/R\n");
    printf("  A: p >= ln(W_MOVE)/ln(k)\n");
    printf("  B: ln(G_MAX)/ln(k) <= p <= ln(1-G_MAX)/ln(1-k)\n\n");

    /* ---- one brush per mesh at the hand-entered exponent ---- */
    printf("p=%.3f metrics\n", P_FIXED);
    printf("%-7s %8s %4s %10s %10s %10s %12s\n",
           "mesh", "h", "n", "affected", "edge_grad", "grad_r", "seam_area");
    for (int d = 0; d < ND; d++) {
        Lattice L;
        if (lat_alloc(&L, hs[d]) != 0) return 1;
        lat_fill(&L, P_FIXED);
        int cells = 0;
        double rg = 0.0;
        double g = edge_gradient(&L, &rg);
        double area = seam_area(&L, &cells);
        printf("%-7s %8.4f %4d %10d %10.4f %10.3f %12.4f\n",
               names[d], L.h, L.n, affected(&L), g, rg, area);
        lat_free(&L);
    }

    /* ---- intervals ---- */
    printf("\nexponent intervals\n");
    printf("%-7s %7s %12s %12s %12s %18s\n",
           "mesh", "k", "A: p_A", "B: p_CORE", "B: p_RIM", "usable");
    double pA_d[ND], pC_d[ND], pR_d[ND];
    for (int d = 0; d < ND; d++) {
        Lattice L;
        if (lat_alloc(&L, hs[d]) != 0) return 1;
        double k = hs[d] / R;
        pA_d[d] = p_coverage(&L);
        pC_d[d] = log(G_MAX) / log(k);
        pR_d[d] = log(1.0 - G_MAX) / log(1.0 - k);
        lat_free(&L);
        double lower = dmax(pA_d[d], pC_d[d]);
        printf("%-7s %7.3f %12.3f %12.3f %12.3f ", names[d], k, pA_d[d], pC_d[d], pR_d[d]);
        if (lower <= pR_d[d]) printf("  [%.3f, %.3f]\n", lower, pR_d[d]);
        else printf("  empty\n");
    }

    /* ---- coarse: the two inequalities disagree ---- */
    printf("\ncoarse conflict\n");
    printf("  A vertex count : p >= %.3f\n", pA_d[0]);
    printf("  B edge gradient: %.3f <= p <= %.3f\n", pC_d[0], pR_d[0]);
    printf("  B is empty on its own (%.3f > %.3f): at small p the centre edge drops\n",
           pC_d[0], pR_d[0]);
    printf("  by k^p > G_MAX, at large p the rim edge drops by 1-(1-k)^p > G_MAX.\n");
    printf("  Intersecting with A: A's floor %.3f is above B's rim ceiling %.3f.\n",
           pA_d[0], pR_d[0]);
    printf("  Sweep p upward: the rim jump exceeds %.2f at p=%.3f, before coverage\n",
           G_MAX, pR_d[0]);
    printf("  reaches %d vertices at p=%.3f. The gradient side breaks first.\n",
           N_MIN, pA_d[0]);
    printf("  Root cause: k=%.2f >= G_MAX=%.2f, so w must fall ~1 -> 0 over only\n",
           hs[0] / R, G_MAX);
    printf("  R/h=%.0f edges. Each such edge carries more than G_MAX for every p.\n", R / hs[0]);

    /* ---- checks: closed forms vs lattice ---- */
    for (int d = 0; d < ND; d++) {
        Lattice L;
        if (lat_alloc(&L, hs[d]) != 0) return 1;
        double k = hs[d] / R;

        /* each boundary must make its own edge sit exactly on G_MAX */
        double gc = jump_core(k, pC_d[d]);
        double gr = jump_rim(k, pR_d[d]);
        if (fabs(gc - G_MAX) > 1e-9) {
            printf("CORE_MISMATCH %s G=%.6f want %.6f\n", names[d], gc, G_MAX);
            failures++;
        }
        if (fabs(gr - G_MAX) > 1e-9) {
            printf("RIM_MISMATCH %s G=%.6f want %.6f\n", names[d], gr, G_MAX);
            failures++;
        }

        /* coverage threshold: on at p_A, off just below */
        lat_fill(&L, pA_d[d]);
        int n_at = affected(&L);
        lat_fill(&L, pA_d[d] * 0.999);
        int n_below = affected(&L);
        if (n_at < N_MIN || n_below >= N_MIN) {
            printf("COVERAGE_MISMATCH %s at=%.4f:%d below=%d\n", names[d], pA_d[d], n_at, n_below);
            failures++;
        }

        /* p_A closed form lands on the first ring, r_cut = k */
        double r_cut = pow(W_MOVE, 1.0 / pA_d[d]);
        if (fabs(r_cut - k) > 2e-3) {
            printf("PCA_MISMATCH %s r_cut=%.4f k=%.4f\n", names[d], r_cut, k);
            failures++;
        }
        lat_free(&L);
    }

    /* ---- placement of the hand-entered p ---- */
    printf("\np=%.3f placement\n", P_FIXED);
    for (int d = 0; d < ND; d++) {
        double lower = dmax(pA_d[d], pC_d[d]);
        int in_usable = (lower <= pR_d[d]) && P_FIXED >= lower && P_FIXED <= pR_d[d];
        int interval_ok = lower <= pR_d[d];
        printf("  %-7s interval %-14s p=%.2f %s\n", names[d],
               interval_ok ? "non-empty" : "empty",
               P_FIXED,
               !interval_ok ? "no usable value" : (in_usable ? "usable" : "too hard, steps"));
    }

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
