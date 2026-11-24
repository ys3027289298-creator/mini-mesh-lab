/* split: connected-component split with per-component index compaction. */
#include "mml_split.h"

#include <string.h>

static int uf_find(int* uf, int x) {
    while (uf[x] != x) {
        uf[x] = uf[uf[x]];
        x = uf[x];
    }
    return x;
}

static void uf_union(int* uf, int a, int b) {
    a = uf_find(uf, a);
    b = uf_find(uf, b);
    if (a != b) uf[b] = a;
}

int mml_split(const MmlMesh* in, MmlComp* out, int maxc) {
    if (!in || !out || maxc <= 0) return -1;
    if (in->nv < 0 || in->nf < 0) return -1;
    if (in->nv > MML_SPLIT_MAXV || in->nf > MML_SPLIT_MAXF) return -1;

    int uf[MML_SPLIT_MAXV], seen[MML_SPLIT_MAXV], g2l[MML_SPLIT_MAXV];

    for (int i = 0; i < in->nv; i++) uf[i] = i;
    for (int i = 0; i < in->nf; i++) {
        const MmlTri* t = &in->f[i];
        if (t->a < 0 || t->b < 0 || t->c < 0 ||
            t->a >= in->nv || t->b >= in->nv || t->c >= in->nv) return -1;
        uf_union(uf, t->a, t->b);
        uf_union(uf, t->b, t->c);
    }

    memset(seen, 0, sizeof seen);
    int nc = 0;
    for (int s = 0; s < in->nv; s++) {
        int r = uf_find(uf, s);
        if (seen[r]) continue;
        seen[r] = 1;
        if (nc >= maxc) return -1;

        MmlComp* c = &out[nc];
        memset(c, 0, sizeof(*c));

        /* compact this component's vertices into its own zero-based array */
        for (int v = 0; v < in->nv; v++) {
            if (uf_find(uf, v) != r) continue;
            g2l[v] = c->m.nv;
            c->src[c->m.nv] = v;
            c->m.v[c->m.nv] = in->v[v];
            c->m.nv++;
        }

        /* renumber this component's faces against its own array only */
        for (int fi = 0; fi < in->nf; fi++) {
            const MmlTri* t = &in->f[fi];
            if (uf_find(uf, t->a) != r) continue;
            c->m.f[c->m.nf].a = g2l[t->a];
            c->m.f[c->m.nf].b = g2l[t->b];
            c->m.f[c->m.nf].c = g2l[t->c];
            c->m.nf++;
        }
        nc++;
    }
    return nc;
}

int mml_check_in_range(const MmlComp* c, int nc) {
    if (!c || nc < 0) return 0;
    for (int k = 0; k < nc; k++)
        for (int i = 0; i < c[k].m.nf; i++) {
            const MmlTri* t = &c[k].m.f[i];
            if (t->a < 0 || t->b < 0 || t->c < 0) return 0;
            if (t->a >= c[k].m.nv || t->b >= c[k].m.nv || t->c >= c[k].m.nv) return 0;
        }
    return 1;
}

int mml_check_dense(const MmlComp* c, int nc) {
    if (!c || nc < 0) return 0;
    if (!mml_check_in_range(c, nc)) return 0;
    for (int k = 0; k < nc; k++) {
        int used[MML_SPLIT_MAXV] = {0};
        for (int i = 0; i < c[k].m.nf; i++) {
            used[c[k].m.f[i].a] = 1;
            used[c[k].m.f[i].b] = 1;
            used[c[k].m.f[i].c] = 1;
        }
        for (int i = 0; i < c[k].m.nv; i++) if (!used[i]) return 0;
    }
    return 1;
}

int mml_check_disjoint(const MmlComp* c, int nc) {
    if (!c || nc < 0) return 0;
    int owner[MML_SPLIT_MAXV];
    for (int i = 0; i < MML_SPLIT_MAXV; i++) owner[i] = -1;
    for (int k = 0; k < nc; k++)
        for (int i = 0; i < c[k].m.nv; i++) {
            int s = c[k].src[i];
            if (s < 0 || s >= MML_SPLIT_MAXV) return 0;
            if (owner[s] != -1) return 0;
            owner[s] = k;
        }
    return 1;
}

int mml_check_conserved(const MmlMesh* in, const MmlComp* c, int nc) {
    if (!in || !c || nc < 0) return 0;
    int nv = 0, nf = 0;
    for (int k = 0; k < nc; k++) {
        nv += c[k].m.nv;
        nf += c[k].m.nf;
    }
    return nv == in->nv && nf == in->nf;
}

/* Translate every local face back through src and match it against the
 * input, one-to-one, preserving winding. */
int mml_check_faces(const MmlMesh* in, const MmlComp* c, int nc) {
    if (!in || !c || nc < 0) return 0;
    if (!mml_check_in_range(c, nc)) return 0;
    int hit[MML_SPLIT_MAXF];
    memset(hit, 0, sizeof hit);
    for (int k = 0; k < nc; k++)
        for (int i = 0; i < c[k].m.nf; i++) {
            const MmlTri* t = &c[k].m.f[i];
            int a = c[k].src[t->a], b = c[k].src[t->b], d = c[k].src[t->c];
            int found = 0;
            for (int j = 0; j < in->nf; j++) {
                if (hit[j]) continue;
                if (in->f[j].a == a && in->f[j].b == b && in->f[j].c == d) {
                    hit[j] = 1;
                    found = 1;
                    break;
                }
            }
            if (!found) return 0;
        }
    for (int j = 0; j < in->nf; j++) if (!hit[j]) return 0;
    return 1;
}
