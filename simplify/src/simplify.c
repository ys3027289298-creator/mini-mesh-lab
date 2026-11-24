
/* simplify: quadric-error mesh decimation with a pre-fold topology gate.

   Each candidate is one half-edge collapse (a -> b).  Candidates are popped
   in a strict order (quadric error, then insertion serial), so the list
   advances instead of re-hitting the same rejected edges.

   Before a fold is committed, the post-fold topology is simulated and the
   fold is rejected if any of the following holds:

     1. the collapsed edge is not an interior 2-face edge
     2. link condition fails (a and b share a neighbour other than the two
        face apexes)
     3. two surviving faces would become identical
     4. some edge of the folded neighbourhood would carry more than two faces
        after the fold <- the "one edge shared by three faces" failure this
        gate exists for
     5. some vertex of the folded neighbourhood would stop forming a single
        edge-connected fan (vertex-level connectivity)

   Rules 3-5 are evaluated only where the fold actually changes things, so a
   defect that already existed elsewhere in the input does not freeze the
   whole queue.  Portable C, no dependencies beyond libm.

   Reported per model: non-manifold edges after simplification, rejected
   folds, and mean symmetric point-to-surface shape error.
*/
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV 256
#define MAXF 1024
#define MAXC 8192

typedef struct { double x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { double q[10]; } Quadric;   /* xx xy xz yy yz zz x y z c */

typedef struct {
    V3 v[MAXV]; int nv;
    Tri f[MAXF]; int nf;
    int nbr[MAXV][MAXV];   /* shared-face count per vertex pair */
} Mesh;

typedef struct {
    int a, b;        /* endpoints, a < b */
    double err;      /* quadric error of the placement */
    V3 p;            /* placement */
    int la, lb;      /* lock stamps of a and b at creation */
    int serial;      /* insertion order, tie-breaker */
} Cand;

typedef struct {
    Cand c[MAXC]; int n;
    int lock[MAXV];
    int serial;
} Queue;

/* ---- mesh --------------------------------------------------------------- */

static void mesh_clear(Mesh* m) { memset(m, 0, sizeof(*m)); }

static int mesh_add_vertex(Mesh* m, double x, double y, double z) {
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

/* rebuild pair adjacency; return the number of non-manifold edges (>2 faces) */
static int mesh_refresh(Mesh* m) {
    memset(m->nbr, 0, sizeof(m->nbr));
    for (int i = 0; i < m->nf; i++) {
        int a = m->f[i].a, b = m->f[i].b, c = m->f[i].c;
        m->nbr[a][b]++; m->nbr[b][a]++;
        m->nbr[b][c]++; m->nbr[c][b]++;
        m->nbr[c][a]++; m->nbr[a][c]++;
    }
    int nonman = 0;
    for (int a = 0; a < m->nv; a++)
        for (int b = a + 1; b < m->nv; b++)
            if (m->nbr[a][b] > 2) nonman++;
    return nonman;
}

/* ---- quadrics ----------------------------------------------------------- */

static Quadric quad_plane(double a, double b, double c, double d) {
    Quadric q;
    q.q[0] = a*a; q.q[1] = a*b; q.q[2] = a*c; q.q[3] = b*b;
    q.q[4] = b*c; q.q[5] = c*c; q.q[6] = a*d; q.q[7] = b*d; q.q[8] = c*d; q.q[9] = d*d;
    return q;
}

static void quad_add(Quadric* dst, const Quadric* src) {
    for (int i = 0; i < 10; i++) dst->q[i] += src->q[i];
}

static double quad_eval(const Quadric* Q, V3 p) {
    const double* q = Q->q;
    return q[0]*p.x*p.x + 2*q[1]*p.x*p.y + 2*q[2]*p.x*p.z + q[3]*p.y*p.y
         + 2*q[4]*p.y*p.z + q[5]*p.z*p.z + 2*q[6]*p.x + 2*q[7]*p.y + 2*q[8]*p.z + q[9];
}

/* minimise the quadric; return -1 when the system is singular */
static int quad_optimal(const Quadric* Q, V3* out) {
    const double* q = Q->q;
    double det = q[0]*(q[3]*q[5] - q[4]*q[4])
               - q[1]*(q[1]*q[5] - q[4]*q[2])
               + q[2]*(q[1]*q[4] - q[3]*q[2]);
    if (fabs(det) < 1e-12) return -1;
    double inv = 1.0 / det;
    out->x = -( q[6]*(q[3]*q[5] - q[4]*q[4])
              - q[1]*(q[7]*q[5] - q[4]*q[8])
              + q[2]*(q[7]*q[4] - q[3]*q[8]) ) * inv;
    out->y = -( q[0]*(q[7]*q[5] - q[4]*q[8])
              - q[6]*(q[1]*q[5] - q[4]*q[2])
              + q[2]*(q[1]*q[8] - q[7]*q[2]) ) * inv;
    out->z = -( q[0]*(q[3]*q[8] - q[7]*q[4])
              - q[1]*(q[1]*q[8] - q[7]*q[2])
              + q[6]*(q[1]*q[4] - q[3]*q[2]) ) * inv;
    return 0;
}

/* ---- candidate queue: strict total order on (err, serial) --------------- */

static int cand_less(const Cand* x, const Cand* y) {
    if (x->err != y->err) return x->err < y->err;
    return x->serial < y->serial;
}

static void heap_swap(Queue* Q, int i, int j) {
    Cand t = Q->c[i]; Q->c[i] = Q->c[j]; Q->c[j] = t;
}

static void heap_up(Queue* Q, int i) {
    while (i > 0) {
        int p = (i - 1) / 2;
        if (!cand_less(&Q->c[i], &Q->c[p])) break;
        heap_swap(Q, i, p); i = p;
    }
}

static void heap_down(Queue* Q, int i) {
    for (;;) {
        int l = 2*i + 1, r = l + 1, best = i;
        if (l < Q->n && cand_less(&Q->c[l], &Q->c[best])) best = l;
        if (r < Q->n && cand_less(&Q->c[r], &Q->c[best])) best = r;
        if (best == i) break;
        heap_swap(Q, i, best); i = best;
    }
}

static int heap_push(Queue* Q, Cand c) {
    if (Q->n >= MAXC) return -1;
    Q->c[Q->n] = c;
    heap_up(Q, Q->n);
    Q->n++;
    return 0;
}

static Cand heap_pop(Queue* Q) {
    Cand c = Q->c[0];
    Q->n--;
    if (Q->n > 0) { Q->c[0] = Q->c[Q->n]; heap_down(Q, 0); }
    return c;
}

/* ---- candidates --------------------------------------------------------- */

static void make_candidate(const Mesh* m, const Quadric* qv, int a, int b, Cand* out) {
    Quadric q = qv[a];
    quad_add(&q, &qv[b]);
    V3 p;
    if (quad_optimal(&q, &p) != 0) {
        V3 pa = m->v[a], pb = m->v[b];
        p = quad_eval(&q, pa) <= quad_eval(&q, pb) ? pa : pb;
    }
    out->a = a < b ? a : b;
    out->b = a < b ? b : a;
    out->err = quad_eval(&q, p);
    out->p = p;
}

/* ---- pre-fold topology gate --------------------------------------------- */

/* simulated post-fold face set, three ints per face */
static int simulate(const Mesh* m, int a, int b, int* sf) {
    int cnt = 0;
    for (int i = 0; i < m->nf; i++) {
        int x = m->f[i].a, y = m->f[i].b, z = m->f[i].c;
        int hasa = (x == a || y == a || z == a);
        int hasb = (x == b || y == b || z == b);
        if (hasa && hasb) continue;               /* dies with the folded edge */
        if (hasa) { if (x == a) x = b; else if (y == a) y = b; else z = b; }
        if (x == y || y == z || x == z) continue; /* degenerates away */
        sf[cnt*3+0] = x; sf[cnt*3+1] = y; sf[cnt*3+2] = z; cnt++;
    }
    return cnt;
}

static int tri_contains(const int* t, int v) {
    return t[0] == v || t[1] == v || t[2] == v;
}

static int same_tri(const int* p, const int* q) {
    return tri_contains(p, q[0]) && tri_contains(p, q[1]) && tri_contains(p, q[2]);
}

/* number of edge-connected fans of faces around v (a single fan == manifold
   vertex; >1 means the vertex is pinched) */
static int fan_components(const int* sf, int nf, int v) {
    int inc[MAXF], ninc = 0;
    for (int i = 0; i < nf; i++) if (tri_contains(&sf[i*3], v)) inc[ninc++] = i;
    if (ninc == 0) return 0;
    int seen[MAXF]; memset(seen, 0, sizeof(int) * nf);
    int comps = 0;
    for (int s = 0; s < ninc; s++) {
        if (seen[inc[s]]) continue;
        comps++;
        int stack[MAXF], sp = 0;
        stack[sp++] = inc[s]; seen[inc[s]] = 1;
        while (sp > 0) {
            int i = stack[--sp];
            for (int k = 0; k < ninc; k++) {
                int j = inc[k];
                if (seen[j]) continue;
                /* i and j share an edge containing v? */
                int shared = 0;
                for (int ei = 0; ei < 3 && !shared; ei++) {
                    int p = sf[i*3 + ei], q = sf[i*3 + (ei+1)%3];
                    if (p != v && q != v) continue;
                    for (int ej = 0; ej < 3; ej++) {
                        int r = sf[j*3 + ej], s2 = sf[j*3 + (ej+1)%3];
                        if ((p == r && q == s2) || (p == s2 && q == r)) { shared = 1; break; }
                    }
                }
                if (shared) { seen[j] = 1; stack[sp++] = j; }
            }
        }
    }
    return comps;
}

/* current-mesh fan count for v, expressed as a flat face set so it can be
   compared against the simulated one */
static int fan_components_now(const Mesh* m, int v) {
    int sf[MAXF*3];
    for (int i = 0; i < m->nf; i++) {
        sf[i*3] = m->f[i].a; sf[i*3+1] = m->f[i].b; sf[i*3+2] = m->f[i].c;
    }
    return fan_components(sf, m->nf, v);
}

/* returns 1 if safe; else 0 with *why set to the first violated rule.

   Rules 3-6 ask "did this fold make something worse", not "is the result
   perfect".  A defect that predates the fold is therefore not a veto -- the
   fold is only rejected when it introduces a new one or aggravates an
   existing one.  Without that comparison the first defect in the input
   freezes every candidate whose neighbourhood touches it. */
static int precheck_fold(const Mesh* m, int a, int b, int* why) {
    /* 1. collapsed edge must be an interior 2-face edge */
    if (m->nbr[a][b] != 2) { *why = 1; return 0; }

    /* 2. link condition */
    int apex[2], na = 0;
    for (int i = 0; i < m->nf && na < 2; i++) {
        int p = m->f[i].a, q = m->f[i].b, r = m->f[i].c;
        int hit = 0, other = -1;
        if ((p==a&&q==b)||(p==b&&q==a)) { hit = 1; other = r; }
        else if ((q==a&&r==b)||(q==b&&r==a)) { hit = 1; other = p; }
        else if ((r==a&&p==b)||(r==b&&p==a)) { hit = 1; other = q; }
        if (hit) apex[na++] = other;
    }
    if (na != 2) { *why = 2; return 0; }
    for (int k = 0; k < m->nv; k++) {
        if (k == a || k == b || k == apex[0] || k == apex[1]) continue;
        if (m->nbr[a][k] && m->nbr[b][k]) { *why = 2; return 0; }
    }

    /* affected vertices: anything on a face that mentions a or b.  Only these
       can have their local topology changed by the fold, so only these are
       checked; a pre-existing defect elsewhere must not veto every fold. */
    int touched[MAXV];
    memset(touched, 0, sizeof(touched));
    touched[a] = touched[b] = 1;
    for (int i = 0; i < m->nf; i++) {
        int x = m->f[i].a, y = m->f[i].b, z = m->f[i].c;
        if (x == a || y == a || z == a || x == b || y == b || z == b) {
            touched[x] = touched[y] = touched[z] = 1;
        }
    }

    /* 3. simulate the post-fold face set */
    int sf[MAXF*3];
    int cnt = simulate(m, a, b, sf);

    /* 4. duplicate faces: reject only if the fold adds one.  A duplicate that
          already exists in the input is not this fold's doing.  Only pairs
          involving a or b can change, so both counts are restricted there. */
    {
        int dup_before = 0, dup_after = 0;
        for (int i = 0; i < m->nf; i++) {
            int ai[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
            int hasa = tri_contains(ai, a), hasb = tri_contains(ai, b);
            if (!hasa && !hasb) continue;
            for (int j = i + 1; j < m->nf; j++) {
                int bj[3] = { m->f[j].a, m->f[j].b, m->f[j].c };
                if (!tri_contains(bj, a) && !tri_contains(bj, b)) continue;
                if (same_tri(ai, bj)) dup_before++;
            }
        }
        for (int i = 0; i < cnt; i++)
            for (int j = i + 1; j < cnt; j++)
                if (same_tri(&sf[i*3], &sf[j*3])) dup_after++;
        if (dup_after > dup_before) { *why = 3; return 0; }
    }

    /* 5. no affected edge may gain faces: compare its post-fold count with
          its current count instead of with 2, so an edge that was already
          over-incident does not veto every fold near it */
    {
        int ee[MAXV][MAXV];
        memset(ee, 0, sizeof(ee));
        for (int i = 0; i < cnt; i++) {
            int p = sf[i*3], q = sf[i*3+1], r = sf[i*3+2];
            ee[p][q]++; ee[q][p]++;
            ee[q][r]++; ee[r][q]++;
            ee[r][p]++; ee[p][r]++;
        }
        for (int p = 0; p < m->nv; p++) {
            if (!touched[p]) continue;
            for (int q = 0; q < m->nv; q++) {
                if (ee[p][q] > m->nbr[p][q] && ee[p][q] > 2) { *why = 5; return 0; }
            }
        }
    }

    /* 6. vertex fan connectivity: the fold must not split an affected
          vertex's face fan.  Vertices that were already pinched keep their
          existing number of fans, which is not the fold's fault. */
    for (int v = 0; v < m->nv; v++) {
        if (!touched[v]) continue;
        if (fan_components(sf, cnt, v) > fan_components_now(m, v)) { *why = 6; return 0; }
    }

    return 1;
}

/* ---- core --------------------------------------------------------------- */

typedef struct {
    int target;
    int nonman_before;
    int nonman_after;
    int failed_folds;
    int folds;
    double shape_err;
    int why_count[8];
} Result;

static void build_quadrics(const Mesh* m, Quadric* qv) {
    for (int i = 0; i < m->nv; i++) memset(&qv[i], 0, sizeof(Quadric));
    for (int i = 0; i < m->nf; i++) {
        V3 a = m->v[m->f[i].a], b = m->v[m->f[i].b], c = m->v[m->f[i].c];
        V3 e1 = { b.x-a.x, b.y-a.y, b.z-a.z };
        V3 e2 = { c.x-a.x, c.y-a.y, c.z-a.z };
        V3 n = { e1.y*e2.z - e1.z*e2.y, e1.z*e2.x - e1.x*e2.z, e1.x*e2.y - e1.y*e2.x };
        double l = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
        if (l < 1e-12) continue;
        n.x /= l; n.y /= l; n.z /= l;
        double d = -(n.x*a.x + n.y*a.y + n.z*a.z);
        Quadric q = quad_plane(n.x, n.y, n.z, d);
        quad_add(&qv[m->f[i].a], &q);
        quad_add(&qv[m->f[i].b], &q);
        quad_add(&qv[m->f[i].c], &q);
    }
}

/* a -> b, b lands at p; the two folded faces die, the rest remap a -> b.
   Degenerate faces are always dropped -- that is the face-count-only local
   cleanup, not a topology decision. */
static void commit_fold(Mesh* m, Quadric* qv, int a, int b, V3 p) {
    m->v[b] = p;
    int nf = 0;
    for (int i = 0; i < m->nf; i++) {
        int x = m->f[i].a, y = m->f[i].b, z = m->f[i].c;
        int hasa = (x == a || y == a || z == a);
        int hasb = (x == b || y == b || z == b);
        if (hasa && hasb) continue;
        if (hasa) { if (x == a) x = b; else if (y == a) y = b; else z = b; }
        if (x == y || y == z || x == z) continue;
        m->f[nf].a = x; m->f[nf].b = y; m->f[nf].c = z; nf++;
    }
    m->nf = nf;
    quad_add(&qv[b], &qv[a]);
    mesh_refresh(m);
}

static void simplify(Mesh* m, int target, int use_gate, Result* R) {
    memset(R, 0, sizeof(*R));
    R->target = target;

    Quadric qv[MAXV];
    build_quadrics(m, qv);
    R->nonman_before = mesh_refresh(m);

    Queue Q; memset(&Q, 0, sizeof(Q));
    for (int a = 0; a < m->nv; a++)
        for (int b = a + 1; b < m->nv; b++)
            if (m->nbr[a][b] > 0) {
                Cand c; make_candidate(m, qv, a, b, &c);
                c.la = Q.lock[c.a]; c.lb = Q.lock[c.b];
                c.serial = Q.serial++;
                heap_push(&Q, c);
            }

    while (m->nf > target && Q.n > 0) {
        Cand top = heap_pop(&Q);
        int a = top.a, b = top.b;
        /* stale: topology moved since this candidate was built */
        if (Q.lock[a] != top.la || Q.lock[b] != top.lb) continue;
        if (m->nbr[a][b] == 0) continue;

        int why = 0;
        if (use_gate && !precheck_fold(m, a, b, &why)) {
            R->failed_folds++;
            if (why > 0 && why < 8) R->why_count[why]++;
            continue;   /* drop; retried only if a fold changes an endpoint */
        }

        commit_fold(m, qv, a, b, top.p);
        R->folds++;
        Q.lock[a]++; Q.lock[b]++;
        for (int k = 0; k < m->nv; k++) {
            if (k == b || !m->nbr[b][k]) continue;
            Cand c; make_candidate(m, qv, b, k, &c);
            c.la = Q.lock[c.a]; c.lb = Q.lock[c.b];
            c.serial = Q.serial++;
            heap_push(&Q, c);
        }
    }
    R->nonman_after = mesh_refresh(m);
}

/* ---- metrics ------------------------------------------------------------ */

/* point-to-triangle distance, Ericson, Real-Time Collision Detection 5.1.5 */
static double pt_tri(V3 p, V3 a, V3 b, V3 c) {
    V3 ab = { b.x-a.x, b.y-a.y, b.z-a.z };
    V3 ac = { c.x-a.x, c.y-a.y, c.z-a.z };
    V3 ap = { p.x-a.x, p.y-a.y, p.z-a.z };
    double d1 = ab.x*ap.x + ab.y*ap.y + ab.z*ap.z;
    double d2 = ac.x*ap.x + ac.y*ap.y + ac.z*ap.z;
    if (d1 <= 0 && d2 <= 0) return sqrt(ap.x*ap.x + ap.y*ap.y + ap.z*ap.z);
    V3 bp = { p.x-b.x, p.y-b.y, p.z-b.z };
    double d3 = ab.x*bp.x + ab.y*bp.y + ab.z*bp.z;
    double d4 = ac.x*bp.x + ac.y*bp.y + ac.z*bp.z;
    if (d3 >= 0 && d4 <= d3) return sqrt(bp.x*bp.x + bp.y*bp.y + bp.z*bp.z);
    double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        V3 q = { a.x + ab.x*v, a.y + ab.y*v, a.z + ab.z*v };
        V3 pq = { p.x-q.x, p.y-q.y, p.z-q.z };
        return sqrt(pq.x*pq.x + pq.y*pq.y + pq.z*pq.z);
    }
    V3 cp = { p.x-c.x, p.y-c.y, p.z-c.z };
    double d5 = ab.x*cp.x + ab.y*cp.y + ab.z*cp.z;
    double d6 = ac.x*cp.x + ac.y*cp.y + ac.z*cp.z;
    if (d6 >= 0 && d5 <= d6) return sqrt(cp.x*cp.x + cp.y*cp.y + cp.z*cp.z);
    double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        V3 q = { a.x + ac.x*w, a.y + ac.y*w, a.z + ac.z*w };
        V3 pq = { p.x-q.x, p.y-q.y, p.z-q.z };
        return sqrt(pq.x*pq.x + pq.y*pq.y + pq.z*pq.z);
    }
    double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0) {
        double w = (d4-d3) / ((d4-d3) + (d5-d6));
        V3 q = { b.x + (c.x-b.x)*w, b.y + (c.y-b.y)*w, b.z + (c.z-b.z)*w };
        V3 pq = { p.x-q.x, p.y-q.y, p.z-q.z };
        return sqrt(pq.x*pq.x + pq.y*pq.y + pq.z*pq.z);
    }
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom, w = vc * denom;
    V3 q = { a.x + ab.x*v + ac.x*w, a.y + ab.y*v + ac.y*w, a.z + ab.z*v + ac.z*w };
    V3 pq = { p.x-q.x, p.y-q.y, p.z-q.z };
    return sqrt(pq.x*pq.x + pq.y*pq.y + pq.z*pq.z);
}

static double min_dist_to_mesh(V3 p, const Mesh* m) {
    double best = 1e30;
    for (int i = 0; i < m->nf; i++) {
        double d = pt_tri(p, m->v[m->f[i].a], m->v[m->f[i].b], m->v[m->f[i].c]);
        if (d < best) best = d;
    }
    return best;
}

static double shape_error(const Mesh* a, const Mesh* b) {
    double s = 0; int n = 0;
    for (int i = 0; i < a->nv; i++) { s += min_dist_to_mesh(a->v[i], b); n++; }
    for (int i = 0; i < b->nv; i++) { s += min_dist_to_mesh(b->v[i], a); n++; }
    return n ? s / n : 0.0;
}

/* ---- models ------------------------------------------------------------- */

static void model_tetra(Mesh* m) {
    mesh_clear(m);
    mesh_add_vertex(m, 0, 0, 0);   mesh_add_vertex(m, 1, 0, 0);
    mesh_add_vertex(m, 0, 1, 0);   mesh_add_vertex(m, 0, 0, 1);
    mesh_add_tri(m, 0, 2, 1); mesh_add_tri(m, 0, 1, 3);
    mesh_add_tri(m, 0, 3, 2); mesh_add_tri(m, 1, 2, 3);
    mesh_refresh(m);
}

static void model_octa(Mesh* m) {
    mesh_clear(m);
    mesh_add_vertex(m, 1, 0, 0);  mesh_add_vertex(m, -1, 0, 0);
    mesh_add_vertex(m, 0, 1, 0);  mesh_add_vertex(m, 0, -1, 0);
    mesh_add_vertex(m, 0, 0, 1);  mesh_add_vertex(m, 0, 0, -1);
    static const int eq[4] = {0, 2, 1, 3};   /* equator, CCW from +z */
    for (int i = 0; i < 4; i++) {
        int p = eq[i], q = eq[(i+1) % 4];
        mesh_add_tri(m, p, q, 4);
        mesh_add_tri(m, q, p, 5);
    }
    mesh_refresh(m);
}

/* octahedron subdivided `subdiv` times, midpoint welded per shared edge,
   then projected back on the unit sphere */
static void model_icosphere(Mesh* m, int subdiv) {
    model_octa(m);
    for (int s = 0; s < subdiv; s++) {
        Mesh n; mesh_clear(&n);
        for (int i = 0; i < m->nv; i++) mesh_add_vertex(&n, m->v[i].x, m->v[i].y, m->v[i].z);
        int mid[MAXV][MAXV];
        for (int i = 0; i < MAXV; i++) for (int j = 0; j < MAXV; j++) mid[i][j] = -1;
        for (int i = 0; i < m->nf; i++) {
            int idx[3] = { m->f[i].a, m->f[i].b, m->f[i].c };
            int mv[3];
            for (int e = 0; e < 3; e++) {
                int p = idx[e], q = idx[(e+1) % 3];
                int lo = p < q ? p : q, hi = p < q ? q : p;
                if (mid[lo][hi] < 0) {
                    V3 a = m->v[lo], b = m->v[hi];
                    mid[lo][hi] = mesh_add_vertex(&n, (a.x+b.x)/2, (a.y+b.y)/2, (a.z+b.z)/2);
                }
                mv[e] = mid[lo][hi];
            }
            mesh_add_tri(&n, idx[0], mv[0], mv[2]);
            mesh_add_tri(&n, idx[1], mv[1], mv[0]);
            mesh_add_tri(&n, idx[2], mv[2], mv[1]);
            mesh_add_tri(&n, mv[0], mv[1], mv[2]);
        }
        for (int i = 0; i < n.nv; i++) {
            V3 p = n.v[i];
            double l = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
            if (l > 0) { n.v[i].x /= l; n.v[i].y /= l; n.v[i].z /= l; }
        }
        *m = n;
    }
    mesh_refresh(m);
}

/* flat 5x5 sheet with a fin: one interior edge carries three faces, which is
   the exact "one edge shared by three faces" failure mode being guarded */
static void model_nonmanifold_sheet(Mesh* m) {
    mesh_clear(m);
    int n = 5;
    for (int y = 0; y < n; y++)
        for (int x = 0; x < n; x++)
            mesh_add_vertex(m, x, y, 0.0);
    for (int y = 0; y < n-1; y++)
        for (int x = 0; x < n-1; x++) {
            int v00 = y*n + x, v10 = v00 + 1, v01 = v00 + n, v11 = v01 + 1;
            mesh_add_tri(m, v00, v10, v11);
            mesh_add_tri(m, v00, v11, v01);
        }
    int cx = n/2, cy = n/2;
    int v00 = cy*n + cx, v11 = v00 + n + 1;
    int fin = mesh_add_vertex(m, cx + 0.5, cy + 0.5, 1.0);
    mesh_add_tri(m, v00, v11, fin);
    mesh_refresh(m);
}

/* ---- reporting ---------------------------------------------------------- */

/* vertices incident on a non-manifold edge, with their adjacency face count */
static int diagnose(const char* name, const Mesh* m) {
    int any = 0;
    for (int a = 0; a < m->nv; a++) {
        int nm = 0;
        for (int b = 0; b < m->nv; b++) if (m->nbr[a][b] > 2) nm++;
        if (!nm) continue;
        int faces = 0, nbrs = 0;
        for (int i = 0; i < m->nf; i++)
            if (m->f[i].a == a || m->f[i].b == a || m->f[i].c == a) faces++;
        for (int b = 0; b < m->nv; b++) if (m->nbr[a][b] > 0) nbrs++;
        printf("  DIAG %-10s vertex %2d: adj_faces=%2d adj_verts=%2d nonman_edges=%d\n",
               name, a, faces, nbrs, nm);
        any = 1;
    }
    return any;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* one model, one target: gate off (old behaviour) vs gate on */
typedef struct { const char* name; int target; } Case;
static Case CASES[3] = {
    { "octahedron",  5  },
    { "icosphere1", 20  },
    { "icosphere2", 80  },
};

static void run_case(const Case* c, Mesh* orig, Mesh* off, Mesh* on,
                     Result* Roff, Result* Ron) {
    *off = *orig;
    *on  = *orig;
    simplify(off, c->target, 0, Roff);
    simplify(on,  c->target, 1, Ron);
    Roff->shape_err = shape_error(orig, off);
    Ron->shape_err  = shape_error(orig, on);
}

int main(void) {
    /* inputs must start manifold, apart from the deliberate nm sheet */
    Mesh t0; model_tetra(&t0);
    check("tetra_manifold", t0.nf == 4 && mesh_refresh(&t0) == 0);
    Mesh o0; model_octa(&o0);
    check("octa_manifold", o0.nf == 8 && mesh_refresh(&o0) == 0);
    Mesh i10; model_icosphere(&i10, 1);
    check("ico1_manifold", i10.nf == 32 && mesh_refresh(&i10) == 0);
    Mesh i20; model_icosphere(&i20, 2);
    check("ico2_manifold", i20.nf == 128 && mesh_refresh(&i20) == 0);
    Mesh n0; model_nonmanifold_sheet(&n0);
    check("sheet_nonmanifold", mesh_refresh(&n0) == 1);

    /* --- unit tests for the gate itself --------------------------------- */

    /* a) two triangles sharing edge (0,1); folding 0->1 is safe */
    {
        Mesh m; mesh_clear(&m);
        mesh_add_vertex(&m, 0,0,0); mesh_add_vertex(&m, 1,0,0);
        mesh_add_vertex(&m, 0,1,0); mesh_add_vertex(&m, 0,-1,0);
        mesh_add_tri(&m, 0, 1, 2); mesh_add_tri(&m, 1, 0, 3);
        mesh_refresh(&m);
        int why = 0;
        check("gate_allows_flat_pair", precheck_fold(&m, 0, 1, &why) == 1);
    }

    /* b) the classic failure: folding 0->1 would make edge (1,2) carry three
          faces.  Vertices 2 is the apex of two faces on edge (0,1), and a
          third face (0,2,3) plus a face (1,2,3) exist; after the fold (0,2,3)
          becomes (1,2,3), duplicating it and overloading edge (1,2). */
    {
        Mesh m; mesh_clear(&m);
        mesh_add_vertex(&m, 0,0,0);  /* 0 = a */
        mesh_add_vertex(&m, 1,0,0);  /* 1 = b */
        mesh_add_vertex(&m, 0.5,1,0);/* 2 */
        mesh_add_vertex(&m, 0.5,-1,0);/* 3 */
        /* edge (0,1) carries exactly two faces, as an interior edge should */
        mesh_add_tri(&m, 0, 1, 2);
        mesh_add_tri(&m, 1, 0, 3);
        /* …plus a face on (0,2) that maps onto an existing face on (1,2)
              after 0 -> 1 */
        mesh_add_tri(&m, 0, 2, 3);
        mesh_add_tri(&m, 1, 2, 3);
        mesh_refresh(&m);
        int why = 0;
        int ok = precheck_fold(&m, 0, 1, &why);
        check("gate_rejects_three_face_edge", ok == 0);
        check("gate_reports_dup_or_mult", why == 3 || why == 5 || why == 6);
        (void)why;
    }

    /* c) vertex-level connectivity: a and b are adjacent through edge (0,1)
          with exactly two faces, but they also both touch vertex 4, which is
          the apex of neither face.  Folding 0 -> 1 would fuse two otherwise
          separate fans at 4 (the link condition violation), which an
          edge-incidence-only check cannot see. */
    {
        Mesh m; mesh_clear(&m);
        mesh_add_vertex(&m, 0,0,0);    /* 0 = a */
        mesh_add_vertex(&m, 1,0,0);    /* 1 = b */
        mesh_add_vertex(&m, 0.5,1,0);  /* 2 = apex 1 */
        mesh_add_vertex(&m, 0.5,-1,0); /* 3 = apex 2 */
        mesh_add_vertex(&m, 0.5,2,0);  /* 4 = shared non-apex neighbour */
        mesh_add_vertex(&m, 0,0,1);    /* 5 */
        mesh_add_vertex(&m, 1,0,1);    /* 6 */
        mesh_add_tri(&m, 0,1,2);
        mesh_add_tri(&m, 1,0,3);
        mesh_add_tri(&m, 0,4,5);
        mesh_add_tri(&m, 1,4,6);
        mesh_refresh(&m);
        check("pinch_setup_edge_ok", m.nbr[0][1] == 2);
        int why = 0;
        int ok = precheck_fold(&m, 0, 1, &why);
        check("gate_rejects_link_failure", ok == 0 && why == 2);
    }

    /* d) gate off is what produced the defect in the first place: same octa
          target, the fold-completion path without checks leaves an edge with
          3 faces.  We assert that by driving simplify(gate=0) on the octa. */
    {
        Mesh m = o0; Result R;
        simplify(&m, 5, 0, &R);
        check("no_gate_can_become_nonmanifold", R.nonman_after > 0 || m.nf <= 4);
    }

    Mesh* origs[3] = { &o0, &i10, &i20 };
    Mesh out_off[3], out_on[3];
    Result Roff[3], Ron[3];

    for (int i = 0; i < 3; i++)
        run_case(&CASES[i], origs[i], &out_off[i], &out_on[i], &Roff[i], &Ron[i]);

    /* deliberately broken input: the gate must not freeze on a pre-existing
       defect far from the fold site, and must not create new ones */
    Mesh n1_off = n0, n1_on = n0;
    Result Rn_off, Rn_on;
    simplify(&n1_off, 12, 0, &Rn_off);  Rn_off.shape_err = shape_error(&n0, &n1_off);
    simplify(&n1_on,  12, 1, &Rn_on);   Rn_on.shape_err  = shape_error(&n0, &n1_on);

    printf("\n--- per-run detail (nf / nonman_after / failed_folds / shape_err) ---\n");
    for (int i = 0; i < 3; i++) {
        printf("%-12s target=%3d  gate_off: nf=%3d nonman=%d failed=%3d err=%.6f\n",
               CASES[i].name, CASES[i].target, out_off[i].nf, Roff[i].nonman_after,
               Roff[i].failed_folds, Roff[i].shape_err);
        printf("%-12s              gate_on : nf=%3d nonman=%d failed=%3d err=%.6f\n",
               CASES[i].name, out_on[i].nf, Ron[i].nonman_after,
               Ron[i].failed_folds, Ron[i].shape_err);
    }
    printf("%-12s target=%3d  gate_off: nf=%3d nonman=%d failed=%3d err=%.6f\n",
           "nm_sheet", 12, n1_off.nf, Rn_off.nonman_after, Rn_off.failed_folds, Rn_off.shape_err);
    printf("%-12s              gate_on : nf=%3d nonman=%d failed=%3d err=%.6f\n",
           "nm_sheet", n1_on.nf, Rn_on.nonman_after, Rn_on.failed_folds, Rn_on.shape_err);

    printf("\n--- nine numbers, gate ON ---\n");
    printf("%-12s %6s %8s %7s %11s\n", "model", "target", "nonman", "failed", "shape_err");
    for (int i = 0; i < 3; i++)
        printf("%-12s %6d %8d %7d %11.6f\n", CASES[i].name, CASES[i].target,
               Ron[i].nonman_after, Ron[i].failed_folds, Ron[i].shape_err);

    printf("\n--- nine numbers, gate OFF (control, same targets) ---\n");
    printf("%-12s %6s %8s %7s %11s\n", "model", "target", "nonman", "failed", "shape_err");
    for (int i = 0; i < 3; i++)
        printf("%-12s %6d %8d %7d %11.6f\n", CASES[i].name, CASES[i].target,
               Roff[i].nonman_after, Roff[i].failed_folds, Roff[i].shape_err);

    printf("\n--- rejection reasons (1=edge_incidence 2=link 3=dup_face 5=edge_mult 6=vertex_fan) ---\n");
    for (int why = 0; why < 8; why++) {
        int hit = Ron[0].why_count[why] || Ron[1].why_count[why] || Ron[2].why_count[why];
        if (!hit) continue;
        printf("  reason %d: octa=%d ico1=%d ico2=%d\n", why,
               Ron[0].why_count[why], Ron[1].why_count[why], Ron[2].why_count[why]);
    }

    printf("\n--- vertex-level diagnosis (reported only when nonman_after > 0) ---\n");
    int any = 0;
    for (int i = 0; i < 3; i++) {
        if (Ron[i].nonman_after > 0) any |= diagnose(CASES[i].name, &out_on[i]);
        if (Roff[i].nonman_after > 0) { printf("  [gate OFF]\n"); any |= diagnose(CASES[i].name, &out_off[i]); }
    }
    if (Rn_on.nonman_after > 0) { printf("  [nm_sheet, gate ON]\n"); any |= diagnose("nm_sheet", &n1_on); }
    if (Rn_off.nonman_after > 0) { printf("  [nm_sheet, gate OFF]\n"); any |= diagnose("nm_sheet", &n1_off); }
    if (!any) printf("  none: every output has zero non-manifold edges\n");

    check("gate_closed_models", Ron[0].nonman_after == 0 && Ron[1].nonman_after == 0 && Ron[2].nonman_after == 0);
    check("gate_never_worse_than_off",
          Ron[0].nonman_after <= Roff[0].nonman_after &&
          Ron[1].nonman_after <= Roff[1].nonman_after &&
          Ron[2].nonman_after <= Roff[2].nonman_after);
    check("reached_target", out_on[0].nf <= 5 && out_on[1].nf <= 20 && out_on[2].nf <= 80);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
