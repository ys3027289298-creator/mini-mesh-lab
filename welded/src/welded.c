
/* welded: vertex welding by position quantization, with face remap. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 64
#define MAXF 64

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;

typedef struct {
    V3 v[MAXV]; int nv;
    Tri f[MAXF]; int nf;
    int map[MAXV];      /* old index -> welded index */
    int welded_count;
    int dropped_degenerate; /* faces removed after welding (two corners coincide) */
    int dropped_invalid;    /* faces skipped because an index was out of range */
} Mesh;

static void mesh_init(Mesh* m) { memset(m, 0, sizeof(*m)); for (int i = 0; i < MAXV; i++) m->map[i] = -1; }

static int mesh_add_vertex(Mesh* m, float x, float y, float z) {
    if (m->nv >= MAXV) return -1;
    m->v[m->nv].x = x; m->v[m->nv].y = y; m->v[m->nv].z = z;
    return m->nv++;
}

static int mesh_add_tri(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAXF) return -1;
    if (a < 0 || b < 0 || c < 0 || a >= m->nv || b >= m->nv || c >= m->nv) return -1;
    if (a == b || b == c || a == c) return -1;
    m->f[m->nf].a = a; m->f[m->nf].b = b; m->f[m->nf].c = c;
    return m->nf++;
}

/* weld vertices closer than eps, remap faces, drop degenerate/invalid.
 *
 * Class admission uses complete linkage: vertex i joins an existing class j
 * only when it is within eps of EVERY earlier vertex already in j. That gives
 * the soundness guarantee  map[u] == map[v]  =>  |pos(u)-pos(v)| <= eps.
 * Consequently a triangle whose own three edges are all longer than eps can
 * never have two corners welded together and be dropped, in any vertex order.
 * (The old first-seed greedy test only compared against the representative, so
 * two points farther than eps apart could share a class via a third point,
 * order-dependently - the intermittent missing face.) */
static int mesh_weld(Mesh* m, float eps) {
    if (eps <= 0.0f) return -1;
    float eps2 = eps * eps;
    V3 uniq[MAXV]; int nu = 0;
    for (int i = 0; i < m->nv; i++) {
        int found = -1;
        for (int j = 0; j < nu && found < 0; j++) {
            int ok = 1;
            for (int k = 0; k < i; k++) {
                if (m->map[k] != j) continue;
                float dx = m->v[i].x - m->v[k].x;
                float dy = m->v[i].y - m->v[k].y;
                float dz = m->v[i].z - m->v[k].z;
                if (dx*dx + dy*dy + dz*dz > eps2) { ok = 0; break; }
            }
            if (ok) found = j;
        }
        if (found < 0) { uniq[nu] = m->v[i]; found = nu++; }
        m->map[i] = found;
    }

    m->dropped_degenerate = 0;
    m->dropped_invalid = 0;
    int nf = 0;
    for (int i = 0; i < m->nf; i++) {
        int ia = m->f[i].a, ib = m->f[i].b, ic = m->f[i].c;
        /* Range-check before indexing map[]: a mesh loaded from disk or built
         * by hand can carry an out-of-range index; the old code indexed
         * map[ia] blindly here (UBSan/ASan: index out of bounds). */
        if (ia < 0 || ib < 0 || ic < 0 ||
            ia >= m->nv || ib >= m->nv || ic >= m->nv) {
            m->dropped_invalid++;
            continue;
        }
        int a = m->map[ia], b = m->map[ib], c = m->map[ic];
        if (a == b || b == c || a == c) { m->dropped_degenerate++; continue; }
        m->f[nf].a = a; m->f[nf].b = b; m->f[nf].c = c; nf++;
    }
    m->nf = nf;
    memcpy(m->v, uniq, sizeof(V3) * nu);
    m->nv = nu;
    m->welded_count = nu;
    return nu;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* next lexicographic permutation of 4 ints; returns 0 when it wraps. */
static int next_perm4(int* a) {
    int i = 2;
    while (i >= 0 && a[i] >= a[i+1]) i--;
    if (i < 0) { for (int k=0;k<4;k++) a[k]=k; return 0; }
    int j = 3;
    while (a[j] <= a[i]) j--;
    int t=a[i]; a[i]=a[j]; a[j]=t;
    for (int l=i+1, r=3; l<r; l++, r--) { t=a[l]; a[l]=a[r]; a[r]=t; }
    return 1;
}

/* Minimal repros for the "which branch actually loses a face" investigation.
 * Every case is deterministic and runs through the real mesh_weld. */
static void branch_regression_checks(void) {
    /* 1) literal duplicate index at insertion: always rejected */
    {
        Mesh d; mesh_init(&d);
        mesh_add_vertex(&d,0,0,0); mesh_add_vertex(&d,1,0,0); mesh_add_vertex(&d,0,1,0);
        check("br_dup_index_rejected", mesh_add_tri(&d,0,0,2) == -1 && d.nf == 0);
    }

    /* 2) out-of-range index at insertion, including index == nv: rejected */
    {
        Mesh o; mesh_init(&o);
        mesh_add_vertex(&o,0,0,0); mesh_add_vertex(&o,1,0,0); mesh_add_vertex(&o,0,1,0);
        int big = mesh_add_tri(&o,0,1,99);
        int edge = mesh_add_tri(&o,0,1,3);   /* == nv, one past the last valid 2 */
        int valid = mesh_add_tri(&o,0,1,2);
        check("br_oob_add_rejected", big == -1 && edge == -1 && valid == 0 && o.nf == 1);
    }

    /* 3) geometrically degenerate (collinear) but 3 distinct indices:
     * there is no area cull, so it survives welding at a small eps. */
    {
        Mesh g; mesh_init(&g);
        mesh_add_vertex(&g,0,0,0); mesh_add_vertex(&g,1,0,0); mesh_add_vertex(&g,2,0,0);
        mesh_add_tri(&g,0,1,2);
        int nu = mesh_weld(&g,1e-4f);
        check("br_geodegen_distinct_survives", nu == 3 && g.nf == 1);
    }

    /* 4) winding conflict: no orientation/cull branch exists, both survive */
    {
        Mesh w; mesh_init(&w);
        mesh_add_vertex(&w,0,0,0); mesh_add_vertex(&w,1,0,0); mesh_add_vertex(&w,0,1,0);
        mesh_add_tri(&w,0,1,2);   /* CCW */
        mesh_add_tri(&w,0,2,1);   /* CW  */
        mesh_weld(&w,1e-4f);
        check("br_winding_both_kept", w.nf == 2 &&
              w.f[0].a==0 && w.f[0].b==1 && w.f[0].c==2 &&
              w.f[1].a==0 && w.f[1].b==2 && w.f[1].c==1);
    }

    /* 5) Regression for the intermittent missing face.
     * Geometry: A=(0,0,0) B=(+.09,0,0) C=(-.09,0,0) R=(0,1,0), eps=.1.
     * |A-B|=|A-C|=.09<=eps but |B-C|=.18>eps, and every edge of BCR is >eps.
     * HEAD's first-seed greedy grouped B,C through A in 8/24 vertex orders and
     * dropped BCR; complete-linkage guarantees equal map ids imply a distance
     * <=eps, so B,C can never share a class here and BCR survives every order. */
    {
        const float px[4]={0.f, .09f,-.09f,0.f};
        const float py[4]={0.f, 0.f, 0.f, 1.f};
        int p[4]={0,1,2,3};
        int lost=0, kept=0, perms=0;
        do {
            Mesh q; mesh_init(&q);
            for (int k=0;k<4;k++) { int L=p[k]; mesh_add_vertex(&q,px[L],py[L],0); }
            int idx[4]; for (int k=0;k<4;k++) idx[p[k]]=k; /* label -> new index */
            mesh_add_tri(&q, idx[1], idx[2], idx[3]);       /* B,C,R */
            mesh_weld(&q,0.1f);
            if (q.nf == 0) lost++; else kept++;
            perms++;
        } while (next_perm4(p));
        printf("    [order-fix] %d permutations: %d drop BCR, %d keep BCR\n",
               perms, lost, kept);
        check("br_order_permutations_count", perms == 24);
        check("br_face_kept_all_orders", lost == 0 && kept == 24);

        Mesh a; mesh_init(&a);   /* A,B,C,R: HEAD dropped this one */
        mesh_add_vertex(&a,0,0,0); mesh_add_vertex(&a,.09f,0,0);
        mesh_add_vertex(&a,-.09f,0,0); mesh_add_vertex(&a,0,1,0);
        mesh_add_tri(&a,1,2,3); mesh_weld(&a,.1f);
        check("br_order_abcr_now_kept", a.nf == 1 && a.dropped_degenerate == 0);

        Mesh b; mesh_init(&b);   /* B,C,A,R: already kept by HEAD */
        mesh_add_vertex(&b,.09f,0,0); mesh_add_vertex(&b,-.09f,0,0);
        mesh_add_vertex(&b,0,0,0); mesh_add_vertex(&b,0,1,0);
        mesh_add_tri(&b,0,1,3); mesh_weld(&b,.1f);
        check("br_order_bcar_still_kept", b.nf == 1);
    }

    /* 6) out-of-range index reaching mesh_weld directly (imported mesh):
     * must be counted and skipped, never index map[] out of bounds. */
    {
        Mesh x; mesh_init(&x);
        mesh_add_vertex(&x,0,0,0); mesh_add_vertex(&x,1,0,0); mesh_add_vertex(&x,0,1,0);
        x.f[0].a=0; x.f[0].b=1; x.f[0].c=99; x.nf=1;
        int nu = mesh_weld(&x,1e-3f);
        check("br_weld_oob_safe", nu == 3 && x.nf == 0 && x.dropped_invalid == 1);
    }

    /* 7) drop accounting: genuinely collapsed face counted, keeper kept */
    {
        Mesh x; mesh_init(&x);
        mesh_add_vertex(&x,0,0,0); mesh_add_vertex(&x,1,0,0);
        mesh_add_vertex(&x,0,1,0); mesh_add_vertex(&x,0,0,0); /* 3 == 0 */
        mesh_add_tri(&x,0,1,2);   /* keeper  */
        mesh_add_tri(&x,0,3,1);   /* collapses */
        mesh_weld(&x,1e-3f);
        check("br_drop_accounting", x.nf == 1 && x.dropped_degenerate == 1 &&
              x.dropped_invalid == 0);
    }

    /* 8) face capacity boundary: MAXF stored, next rejected (no off-by-one) */
    {
        Mesh c; mesh_init(&c);
        mesh_add_vertex(&c,0,0,0); mesh_add_vertex(&c,1,0,0); mesh_add_vertex(&c,0,1,0);
        int accepted=0;
        for (int i=0;i<MAXF+1;i++) if (mesh_add_tri(&c,0,1,2)==i) accepted++;
        check("br_face_capacity", accepted == MAXF && c.nf == MAXF);
    }
}

int main(void) {
    Mesh m;
    mesh_init(&m);
    mesh_add_vertex(&m, 0, 0, 0);
    mesh_add_vertex(&m, 1, 0, 0);
    mesh_add_vertex(&m, 0, 1, 0);
    mesh_add_vertex(&m, 0, 0, 0);        /* duplicate of 0 */
    mesh_add_vertex(&m, 1.0001f, 0, 0);  /* near duplicate of 1 */
    check("added", m.nv == 5);

    check("tri_ok", mesh_add_tri(&m, 0, 1, 2) == 0);
    check("tri_dup_rejected", mesh_add_tri(&m, 0, 0, 2) == -1);
    check("tri_range_rejected", mesh_add_tri(&m, 0, 1, 99) == -1);
    mesh_add_tri(&m, 3, 4, 2);           /* welds to same as first: survives */
    mesh_add_tri(&m, 0, 3, 1);           /* 0 and 3 weld together: degenerate, dropped */

    int nu = mesh_weld(&m, 1e-3f);
    check("weld_count", nu == 3);
    check("map_zero", m.map[0] == m.map[3]);
    check("map_one", m.map[1] == m.map[4]);
    check("degenerate_dropped", m.nf == 2);

    check("bad_eps", mesh_weld(&m, 0.0f) == -1);

    Mesh e;
    mesh_init(&e);
    mesh_add_vertex(&e, 0,0,0); mesh_add_vertex(&e, 1,0,0); mesh_add_vertex(&e, 0,1,0);
    mesh_add_tri(&e, 0, 1, 2);
    check("no_weld_needed", mesh_weld(&e, 1e-4f) == 3 && e.nf == 1);

    branch_regression_checks();

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
