/* split test: regression for the component index offset, checked at 2, 3
 * and 5 components, each time verifying three properties together:
 *   - every component's indices start at zero and stay inside its array
 *   - no vertex is shared between components
 *   - vertex and face totals are conserved across the split
 * plus a face-level round trip through the local->source table.
 */
#include "mml_split.h"

#include <stdio.h>
#include <string.h>

#define MAXC 16

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void mesh_init(MmlMesh* m) { memset(m, 0, sizeof(*m)); }

static int av(MmlMesh* m, float x, float y, float z) {
    if (m->nv >= MML_SPLIT_MAXV) return -1;
    m->v[m->nv].x = x; m->v[m->nv].y = y; m->v[m->nv].z = z;
    return m->nv++;
}

static int af(MmlMesh* m, int a, int b, int c) {
    if (m->nf >= MML_SPLIT_MAXF) return -1;
    if (a < 0 || b < 0 || c < 0 || a >= m->nv || b >= m->nv || c >= m->nv) return -1;
    if (a == b || b == c || a == c) return -1;
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* One component per entry of sizes[], component k offset to x = 10k and
 * triangulated as a fan around its first vertex. Sizes differ so a running
 * global vertex base would drift differently per component. */
static void build(MmlMesh* m, const int* sizes, int ncomp) {
    mesh_init(m);
    for (int k = 0; k < ncomp; k++) {
        float x = 10.0f * (float)k;
        int base = m->nv;
        for (int i = 0; i < sizes[k]; i++)
            av(m, x + (float)(i % 3), (float)(i / 3), (float)(i == 1));
        for (int i = 1; i + 1 < sizes[k]; i++)
            af(m, base, base + i, base + i + 1);
    }
}

/* Same components, but their source vertices are interleaved: component k
 * owns k, k+ncomp, k+2*ncomp, ... This separates "shifted by global base"
 * from "not reordered at all". */
static void build_interleaved(MmlMesh* m, const int* sizes, int ncomp) {
    mesh_init(m);
    int total = 0, maxsz = 0, slot[MAXC][MML_SPLIT_MAXV];
    for (int k = 0; k < ncomp; k++) {
        total += sizes[k];
        if (sizes[k] > maxsz) maxsz = sizes[k];
    }
    /* round-robin assignment so neighbouring indices belong to different
     * components even when the components have different sizes */
    int next = 0;
    for (int i = 0; i < maxsz; i++)
        for (int k = 0; k < ncomp; k++)
            if (i < sizes[k]) slot[k][i] = next++;
    for (int i = 0; i < total; i++) av(m, (float)i, (float)(i % 3), 0.0f);
    for (int k = 0; k < ncomp; k++)
        for (int i = 1; i + 1 < sizes[k]; i++)
            af(m, slot[k][0], slot[k][i], slot[k][i + 1]);
}

static void verify(const char* tag, const MmlMesh* in, int ncomp, const int* sizes) {
    MmlComp c[MAXC];
    char name[128];
    int nc = mml_split(in, c, MAXC);
    snprintf(name, sizeof name, "%s_components", tag);
    check(name, nc == ncomp);

    int starts_at_zero = 1, all_used = 1;
    for (int k = 0; k < nc && nc > 0; k++) {
        if (c[k].m.nf > 0) {
            /* the fan's first face is built from the component's first
             * vertex, so its local a must be 0 in every component */
            if (c[k].m.f[0].a != 0) starts_at_zero = 0;
        }
        int used[MML_SPLIT_MAXV] = {0};
        for (int i = 0; i < c[k].m.nf; i++)
            used[c[k].m.f[i].a] = used[c[k].m.f[i].b] = used[c[k].m.f[i].c] = 1;
        for (int i = 0; i < c[k].m.nv; i++) if (!used[i]) all_used = 0;
    }
    snprintf(name, sizeof name, "%s_indices_start_at_zero", tag);
    check(name, starts_at_zero);
    snprintf(name, sizeof name, "%s_in_range", tag);
    check(name, mml_check_in_range(c, nc));
    snprintf(name, sizeof name, "%s_dense", tag);
    check(name, all_used && mml_check_dense(c, nc));
    snprintf(name, sizeof name, "%s_disjoint", tag);
    check(name, mml_check_disjoint(c, nc));
    snprintf(name, sizeof name, "%s_conserved", tag);
    check(name, mml_check_conserved(in, c, nc));
    snprintf(name, sizeof name, "%s_faces_roundtrip", tag);
    check(name, mml_check_faces(in, c, nc));

    if (nc == ncomp) {
        int per_comp = 1, sum_nv = 0, sum_nf = 0;
        for (int k = 0; k < nc; k++) {
            if (c[k].m.nv != sizes[k]) per_comp = 0;
            sum_nv += c[k].m.nv;
            sum_nf += c[k].m.nf;
        }
        snprintf(name, sizeof name, "%s_nv_per_comp", tag);
        check(name, per_comp);
        snprintf(name, sizeof name, "%s_totals", tag);
        check(name, sum_nv == in->nv && sum_nf == in->nf);
    }
}

/* Control: re-adding a running global base must trip the same checks, so a
 * passing run above is not vacuous. */
static void verify_bug_is_caught(const MmlMesh* in, const char* tag) {
    MmlComp c[MAXC];
    char name[128];
    int nc = mml_split(in, c, MAXC);
    int base = 0;
    for (int k = 0; k < nc; k++) {
        for (int i = 0; i < c[k].m.nf; i++) {
            c[k].m.f[i].a += base;
            c[k].m.f[i].b += base;
            c[k].m.f[i].c += base;
        }
        base += c[k].m.nv;
    }
    snprintf(name, sizeof name, "%s_bug_in_range", tag);
    check(name, !mml_check_in_range(c, nc));
    snprintf(name, sizeof name, "%s_bug_faces", tag);
    check(name, !mml_check_faces(in, c, nc));
}

/* Control: sharing one vertex table between components must trip disjoint. */
static void verify_shared_is_caught(const MmlMesh* in, const char* tag) {
    MmlComp c[MAXC];
    char name[128];
    int nc = mml_split(in, c, MAXC);
    for (int k = 1; k < nc; k++)
        for (int i = 0; i < c[k].m.nv; i++) c[k].src[i] = c[0].src[i % c[0].m.nv];
    snprintf(name, sizeof name, "%s_shared_disjoint", tag);
    check(name, !mml_check_disjoint(c, nc));
}

int main(void) {
    const int two[2] = {3, 3};
    const int three[3] = {3, 4, 3};
    const int five[5] = {3, 5, 4, 3, 6};

    MmlMesh m;

    /* minimal repro from the report: two triangles */
    mesh_init(&m);
    av(&m, 0,0,0); av(&m, 1,0,0); av(&m, 0,1,0);
    av(&m, 10,0,0); av(&m, 11,0,0); av(&m, 10,1,0);
    af(&m, 0,1,2); af(&m, 3,4,5);
    {
        MmlComp c[MAXC];
        int nc = mml_split(&m, c, MAXC);
        check("minimal_two_comps", nc == 2);
        check("minimal_comp0_faces", c[0].m.f[0].a == 0 && c[0].m.f[0].b == 1 && c[0].m.f[0].c == 2);
        check("minimal_comp1_faces", c[1].m.f[0].a == 0 && c[1].m.f[0].b == 1 && c[1].m.f[0].c == 2);
        check("minimal_comp1_src", c[1].src[0] == 3 && c[1].src[1] == 4 && c[1].src[2] == 5);
        check("minimal_all_checks",
              mml_check_in_range(c, nc) && mml_check_disjoint(c, nc) &&
              mml_check_conserved(&m, c, nc) && mml_check_faces(&m, c, nc));
    }

    build(&m, two, 2);
    verify("two", &m, 2, two);
    build(&m, three, 3);
    verify("three", &m, 3, three);
    build(&m, five, 5);
    verify("five", &m, 5, five);

    build_interleaved(&m, two, 2);
    verify("two_interleaved", &m, 2, two);
    build_interleaved(&m, three, 3);
    verify("three_interleaved", &m, 3, three);
    build_interleaved(&m, five, 5);
    verify("five_interleaved", &m, 5, five);

    verify_bug_is_caught(&m, "five_interleaved");
    verify_shared_is_caught(&m, "five_interleaved");

    /* a single component must come through untouched */
    {
        MmlMesh one;
        mesh_init(&one);
        av(&one, 0,0,0); av(&one, 1,0,0); av(&one, 0,1,0);
        af(&one, 0,1,2);
        MmlComp c[MAXC];
        int nc = mml_split(&one, c, MAXC);
        check("single_component", nc == 1);
        check("single_unchanged",
              c[0].m.nv == 3 && c[0].m.nf == 1 &&
              c[0].m.f[0].a == 0 && c[0].m.f[0].b == 1 && c[0].m.f[0].c == 2);
    }

    /* an isolated vertex with no face still forms its own component */
    {
        MmlMesh loose;
        mesh_init(&loose);
        av(&loose, 0,0,0); av(&loose, 1,0,0); av(&loose, 0,1,0);
        av(&loose, 5,5,5);
        af(&loose, 0,1,2);
        MmlComp c[MAXC];
        int nc = mml_split(&loose, c, MAXC);
        check("isolated_vertex_component", nc == 2);
        check("isolated_vertex_conserved", mml_check_conserved(&loose, c, nc));
    }

    /* edge cases */
    {
        MmlMesh e;
        mesh_init(&e);
        MmlComp c[MAXC];
        check("empty_zero_comps", mml_split(&e, c, MAXC) == 0);

        MmlMesh t;
        mesh_init(&t);
        av(&t, 0,0,0); av(&t, 1,0,0); av(&t, 0,1,0);
        af(&t, 0,1,2);
        av(&t, 10,0,0); av(&t, 11,0,0); av(&t, 10,1,0);
        af(&t, 3,4,5);
        check("too_few_slots", mml_split(&t, c, 1) == -1);

        MmlMesh bad;
        mesh_init(&bad);
        av(&bad, 0,0,0); av(&bad, 1,0,0); av(&bad, 0,1,0);
        bad.f[0].a = 0; bad.f[0].b = 1; bad.f[0].c = 9; bad.nf = 1;
        check("bad_face_index", mml_split(&bad, c, MAXC) == -1);

        check("null_rejected", mml_split(NULL, c, MAXC) == -1);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
