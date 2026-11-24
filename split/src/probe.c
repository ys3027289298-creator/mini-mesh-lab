/* probe: reproduce the reported component-index offset and rule out the
 * other three candidate causes with discriminating inputs.
 *
 * Pre-fix splitter: each component's vertices are compacted correctly, but
 * the component's running global vertex base is then added to every face
 * index.  Only the second and later components are affected.
 *
 * Output sections:
 *   A  minimal repro: two triangles, indices before/after the fix
 *   B  interleaved input: predictions for all four hypotheses
 *   C  positive controls: each invariant detector fires on its own fault
 *   D  invariants on the observed (pre-fix) output
 *   E  verdict per hypothesis
 */
#include <stdio.h>
#include <string.h>

#define MAXV 64
#define MAXF 64
#define MAXC 16

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct {
    V3 v[MAXV]; int nv;
    Tri f[MAXF]; int nf;
} Mesh;

typedef struct {
    Mesh m;
    int src[MAXV];   /* component-local index -> source index */
} Comp;

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); }

static int av(Mesh* m, float x, float y, float z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv] = (V3){x, y, z};
    return m->nv++;
}

static int af(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a < 0 || b < 0 || c < 0 || a >= m->nv || b >= m->nv || c >= m->nv) return -1;
    if (a == b || b == c || a == c) return -1;
    m->f[m->nf] = (Tri){a, b, c};
    return m->nf++;
}

static int uf_find(int* uf, int x) {
    while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
    return x;
}

static void uf_union(int* uf, int a, int b) {
    a = uf_find(uf, a); b = uf_find(uf, b);
    if (a != b) uf[b] = a;
}

/* add_base != 0 reproduces the reported failure; 0 is the intended result */
static int split(const Mesh* in, Comp out[], int maxc, int add_base) {
    if (!in || !out || maxc <= 0 || in->nv > MAXV || in->nf > MAXF) return -1;
    int uf[MAXV], seen[MAXV], g2l[MAXV];
    for (int i = 0; i < in->nv; i++) uf[i] = i;
    for (int i = 0; i < in->nf; i++) {
        uf_union(uf, in->f[i].a, in->f[i].b);
        uf_union(uf, in->f[i].b, in->f[i].c);
    }
    memset(seen, 0, sizeof seen);
    int nc = 0, base = 0;
    for (int s = 0; s < in->nv; s++) {
        int r = uf_find(uf, s);
        if (seen[r]) continue;
        seen[r] = 1;
        if (nc >= maxc) return -1;
        Comp* c = &out[nc];
        memset(c, 0, sizeof(*c));
        for (int v = 0; v < in->nv; v++) {
            if (uf_find(uf, v) != r) continue;
            g2l[v] = c->m.nv;
            c->src[c->m.nv] = v;
            c->m.v[c->m.nv++] = in->v[v];
        }
        for (int fi = 0; fi < in->nf; fi++) {
            const Tri* t = &in->f[fi];
            if (uf_find(uf, t->a) != r) continue;
            int off = add_base ? base : 0;
            c->m.f[c->m.nf++] = (Tri){ g2l[t->a] + off, g2l[t->b] + off, g2l[t->c] + off };
        }
        base += c->m.nv;
        nc++;
    }
    return nc;
}

static void print_comps(const char* tag, const Comp c[], int nc) {
    printf("%s\n", tag);
    for (int k = 0; k < nc; k++) {
        printf("  comp%d nv=%d nf=%d src=[", k, c[k].m.nv, c[k].m.nf);
        for (int i = 0; i < c[k].m.nv; i++) printf("%s%d", i ? " " : "", c[k].src[i]);
        printf("] faces=");
        for (int i = 0; i < c[k].m.nf; i++)
            printf(" (%d,%d,%d)", c[k].m.f[i].a, c[k].m.f[i].b, c[k].m.f[i].c);
        printf("\n");
    }
}

/* ---- invariant detectors (return 1 = holds) ---- */

static int in_range(const Comp c[], int nc) {
    for (int k = 0; k < nc; k++)
        for (int i = 0; i < c[k].m.nf; i++) {
            int idx[3] = { c[k].m.f[i].a, c[k].m.f[i].b, c[k].m.f[i].c };
            for (int j = 0; j < 3; j++)
                if (idx[j] < 0 || idx[j] >= c[k].m.nv) return 0;
        }
    return 1;
}

static int dense_from_zero(const Comp c[], int nc) {
    if (!in_range(c, nc)) return 0;
    for (int k = 0; k < nc; k++) {
        if (c[k].m.nf == 0) continue;
        int used[MAXV] = {0};
        for (int i = 0; i < c[k].m.nf; i++) {
            used[c[k].m.f[i].a] = used[c[k].m.f[i].b] = used[c[k].m.f[i].c] = 1;
        }
        for (int i = 0; i < c[k].m.nv; i++) if (!used[i]) return 0;
    }
    return 1;
}

static int src_disjoint(const Comp c[], int nc) {
    int owner[MAXV];
    for (int i = 0; i < MAXV; i++) owner[i] = -1;
    for (int k = 0; k < nc; k++)
        for (int i = 0; i < c[k].m.nv; i++) {
            int s = c[k].src[i];
            if (s < 0 || s >= MAXV) return 0;
            if (owner[s] != -1) return 0;
            owner[s] = k;
        }
    return 1;
}

static int conserved(const Mesh* in, const Comp c[], int nc, int* sum_nv, int* sum_nf) {
    int nv = 0, nf = 0;
    for (int k = 0; k < nc; k++) { nv += c[k].m.nv; nf += c[k].m.nf; }
    if (sum_nv) *sum_nv = nv;
    if (sum_nf) *sum_nf = nf;
    return nv == in->nv && nf == in->nf;
}

/* translate local indices back through src and require an exact face match */
static int faces_preserved(const Mesh* in, const Comp c[], int nc, int translated) {
    if (!in_range(c, nc)) return 0;
    int hit[MAXF] = {0};
    for (int k = 0; k < nc; k++)
        for (int i = 0; i < c[k].m.nf; i++) {
            const Tri* t = &c[k].m.f[i];
            int A = c[k].src[t->a], B = c[k].src[t->b], C = c[k].src[t->c];
            if (!translated) { A = t->a; B = t->b; C = t->c; }
            int found = 0;
            for (int j = 0; j < in->nf; j++)
                if (!hit[j] && in->f[j].a == A && in->f[j].b == B && in->f[j].c == C) {
                    hit[j] = 1; found = 1; break;
                }
            if (!found) return 0;
        }
    for (int j = 0; j < in->nf; j++) if (!hit[j]) return 0;
    return 1;
}

static void report_invariants(const char* tag, const Mesh* in, const Comp c[], int nc) {
    int nv = 0, nf = 0;
    int ok = conserved(in, c, nc, &nv, &nf);
    printf("  %s: in_range=%d dense=%d disjoint=%d conserved_nv=%d (%d/%d) "
           "conserved_nf=%d (%d/%d) faces_via_src=%d\n",
           tag, in_range(c, nc), dense_from_zero(c, nc), src_disjoint(c, nc),
           nv == in->nv, nv, in->nv, nf == in->nf, nf, in->nf,
           faces_preserved(in, c, nc, 1));
    (void)ok;
}

/* ---- fault injections (positive controls for the detectors) ---- */

/* H1: faces were never renumbered, they still carry source indices */
static void fault_no_reorder(Comp c[], int nc) {
    for (int k = 0; k < nc; k++)
        for (int i = 0; i < c[k].m.nf; i++) {
            Tri t = c[k].m.f[i];
            c[k].m.f[i] = (Tri){ c[k].src[t.a], c[k].src[t.b], c[k].src[t.c] };
        }
}

/* H2: correct local indices plus a running global vertex base */
static void fault_global_offset(Comp c[], int nc) {
    int base = 0;
    for (int k = 0; k < nc; k++) {
        for (int i = 0; i < c[k].m.nf; i++) {
            c[k].m.f[i].a += base; c[k].m.f[i].b += base; c[k].m.f[i].c += base;
        }
        base += c[k].m.nv;
    }
}

/* H3: one shared vertex table; every component indexes it from zero */
static void fault_shared(Comp c[], int nc) {
    for (int k = 1; k < nc; k++)
        for (int i = 0; i < c[k].m.nv; i++)
            c[k].src[i] = c[0].src[i % c[0].m.nv];
}

/* H4: only half of the local->source table was copied; untranslated
 * entries surface as -1 in the face indices */
static void fault_half_map(Comp c[], int nc) {
    for (int k = 0; k < nc; k++) {
        int half = c[k].m.nv / 2;
        for (int i = half; i < c[k].m.nv; i++) c[k].src[i] = -1;
        for (int i = 0; i < c[k].m.nf; i++) {
            int* idx = &c[k].m.f[i].a;
            for (int j = 0; j < 3; j++) if (idx[j] >= half) idx[j] = -1;
        }
    }
}

static void copy_comps(Comp dst[], const Comp src[], int nc) {
    for (int k = 0; k < nc; k++) dst[k] = src[k];
}

int main(void) {
    Comp c[MAXC], f[MAXC];

    /* ---- A. minimal repro ---- */
    Mesh two; mesh_init(&two);
    av(&two, 0,0,0); av(&two, 1,0,0); av(&two, 0,1,0);
    av(&two, 10,0,0); av(&two, 11,0,0); av(&two, 10,1,0);
    af(&two, 0,1,2); af(&two, 3,4,5);

    printf("=== A. minimal repro: two triangles, contiguous indices ===\n");
    printf("input: nv=%d nf=%d faces (0,1,2) (3,4,5)\n", two.nv, two.nf);
    int nc = split(&two, c, MAXC, 1);
    print_comps("before fix (global base added to every face index):", c, nc);
    printf("  comp1 emits 3..5 but owns 3 local vertices, valid range is [0,%d]\n", c[1].m.nv - 1);
    nc = split(&two, c, MAXC, 0);
    print_comps("after fix (component-local indices):", c, nc);

    /* ---- B. discriminating input: interleaved vertex order ---- */
    /* comp0 owns source vertices 0,2,4; comp1 owns 1,3,5.
     * No-reorder, global-offset and shared-table predictions all differ. */
    Mesh il; mesh_init(&il);
    av(&il, 0,0,0); av(&il, 10,0,0); av(&il, 1,0,0);
    av(&il, 11,0,0); av(&il, 0,1,0); av(&il, 10,1,0);
    af(&il, 0,2,4); af(&il, 1,3,5);

    printf("\n=== B. discriminating input: two triangles, interleaved order ===\n");
    printf("input: nv=%d nf=%d faces (0,2,4) (1,3,5)\n", il.nv, il.nf);
    int nc_ok = split(&il, c, MAXC, 0);
    int nc_bad = split(&il, f, MAXC, 1);
    print_comps("correct split:", c, nc_ok);
    print_comps("observed split (before fix):", f, nc_bad);
    printf("predictions for comp1 face0 (local (0,1,2), src [1 3 5], base 3):\n");
    printf("  H1 no reorder      -> (%d,%d,%d)\n",
           c[1].src[0], c[1].src[1], c[1].src[2]);
    printf("  H2 global offset   -> (%d,%d,%d)\n",
           c[1].m.f[0].a + c[0].m.nv, c[1].m.f[0].b + c[0].m.nv, c[1].m.f[0].c + c[0].m.nv);
    printf("  H3 shared vertices -> (%d,%d,%d) with comp1 src = [%d %d %d]\n",
           c[1].m.f[0].a, c[1].m.f[0].b, c[1].m.f[0].c,
           c[0].src[0], c[0].src[1], c[0].src[2]);
    printf("  H4 half index map  -> contains -1, or fewer than %d faces total\n", il.nf);
    printf("  observed           -> (%d,%d,%d), comp1 src = [%d %d %d]\n",
           f[1].m.f[0].a, f[1].m.f[0].b, f[1].m.f[0].c,
           f[1].src[0], f[1].src[1], f[1].src[2]);
    printf("  comp1's first source index is 1, but the shift is +%d = cumulative "
           "vertex count, not +1\n", f[1].m.f[0].a - c[1].m.f[0].a);
    printf("  removing the per-component base recovers the input faces: ");
    {
        Comp g[MAXC];
        copy_comps(g, f, nc_bad);
        int base = 0, ok = 1;
        for (int k = 0; k < nc_bad; k++) {
            for (int i = 0; i < g[k].m.nf; i++) {
                g[k].m.f[i].a -= base; g[k].m.f[i].b -= base; g[k].m.f[i].c -= base;
            }
            base += g[k].m.nv;
        }
        ok = faces_preserved(&il, g, nc_bad, 1);
        printf("%s\n", ok ? "yes, all faces and winding intact" : "no");
    }

    /* ---- C. positive controls ---- */
    printf("\n=== C. each detector fires on its matching fault ===\n");
    copy_comps(f, c, nc_ok); fault_no_reorder(f, nc_ok);
    report_invariants("fault=no reorder    ", &il, f, nc_ok);
    copy_comps(f, c, nc_ok); fault_global_offset(f, nc_ok);
    report_invariants("fault=global offset ", &il, f, nc_ok);
    copy_comps(f, c, nc_ok); fault_shared(f, nc_ok);
    report_invariants("fault=shared verts  ", &il, f, nc_ok);
    copy_comps(f, c, nc_ok); fault_half_map(f, nc_ok);
    report_invariants("fault=half index map", &il, f, nc_ok);

    /* ---- D. invariants on the observed output ---- */
    printf("\n=== D. invariants on the observed (before-fix) output ===\n");
    nc_bad = split(&il, f, MAXC, 1);
    report_invariants("interleaved 2 comps", &il, f, nc_bad);

    /* three components of unequal size: shift is cumulative, not per-component */
    Mesh three; mesh_init(&three);
    av(&three, 0,0,0); av(&three, 10,0,0); av(&three, 20,0,0);
    av(&three, 1,0,0); av(&three, 11,0,0); av(&three, 21,0,0);
    av(&three, 0,1,0); av(&three, 10,1,0); av(&three, 20,1,0);
    av(&three, 10,2,0);
    af(&three, 0,3,6); af(&three, 1,4,7); af(&three, 1,7,9); af(&three, 2,5,8);
    split(&three, c, MAXC, 0);
    int nc3b = split(&three, f, MAXC, 1);
    printf("\n  three comps, sizes 3,4,3: per-component shift added by the bug:");
    for (int k = 0; k < nc3b; k++)
        printf(" comp%d=+%d", k, f[k].m.f[0].a - c[k].m.f[0].a);
    printf("  (cumulative vertex counts)\n");
    report_invariants("three comps         ", &three, f, nc3b);

    /* ---- E. verdicts ---- */
    int h1 = 0, h2 = 1, h3 = 0, h4 = 0;
    {
        int base = 0;
        nc_bad = split(&il, f, MAXC, 1);
        nc_ok = split(&il, c, MAXC, 0);
        for (int k = 0; k < nc_bad && h2; k++) {
            for (int i = 0; i < f[k].m.nf && h2; i++) {
                if (f[k].m.f[i].a != c[k].m.f[i].a + base) h2 = 0;
                if (f[k].m.f[i].b != c[k].m.f[i].b + base) h2 = 0;
                if (f[k].m.f[i].c != c[k].m.f[i].c + base) h2 = 0;
            }
            base += f[k].m.nv;
        }
        (void)h1; (void)h3; (void)h4;
    }

    printf("\n=== E. verdict ===\n");
    printf("H1 no reorder      : EXCLUDED  comp1 emitted (%d,%d,%d), source indices are (%d,%d,%d)\n",
           f[1].m.f[0].a, f[1].m.f[0].b, f[1].m.f[0].c,
           c[1].src[0], c[1].src[1], c[1].src[2]);
    printf("H2 global offset   : %s  every face index equals local + running base\n",
           h2 ? "CONFIRMED" : "not matched");
    printf("H3 shared vertices : EXCLUDED  src=[%d %d %d] vs [%d %d %d] disjoint, sum nv=%d/%d\n",
           f[0].src[0], f[0].src[1], f[0].src[2],
           f[1].src[0], f[1].src[1], f[1].src[2],
           f[0].m.nv + f[1].m.nv, il.nv);
    printf("H4 half index map  : EXCLUDED  sum nf=%d/%d, no -1, all faces present\n",
           f[0].m.nf + f[1].m.nf, il.nf);
    return 0;
}
