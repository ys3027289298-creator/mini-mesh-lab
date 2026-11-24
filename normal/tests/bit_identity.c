
/* bit_identity: prove the U1..U4 refactor is bit-for-bit identical to the
 * pre-split implementation (tests/normal_ref.c, a symbol-renamed copy of
 * git HEAD:normal/src/normal.c).  Both implementations are compiled from
 * separate translation units and called on identical inputs.
 *
 * Comparison policy per component:
 *  - finite components: exact bit comparison (memcmp of the 32-bit pattern)
 *  - NaN components:   NaN-ness only. NaN payload propagation through SSE
 *                       arithmetic is not guaranteed by the C standard and,
 *                       empirically, already differs between two builds of
 *                       the SAME monolithic source because inlining changes
 *                       instruction selection; no source-level refactor can
 *                       promise payload identity. Every finite bit still must
 *                       match, and any NaN-vs-finite disagreement is a fail.
 */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../src/normal.h"

int vertex_normals_ref(const Mesh* m, V3* out);

static uint32_t rng_state = 0x12345678u;
static uint32_t rng_next(void) {
    /* xorshift32 */
    uint32_t x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_state = x;
    return x;
}

/* stress float space: mostly ordinary finite values, with rare signed
 * zeros, infinities, subnormals and NaNs. Keeping the special classes rare
 * means the large meshes (up to 32 verts * 3 components) still usually
 * contain no NaN at all, so the exact-bit path gets most of the traffic.
 * When finite_only is set, NaN/Inf classes are remapped to finite ones. */
static int finite_only = 0;
static float stress_float(void) {
    uint32_t bits = rng_next();
    switch (bits & 0x1Fu) {                 /* /32 distribution */
        case 0: {                           /* subnormal magnitude */
            uint32_t b = (rng_next() & 0x007FFFFFu) | ((bits & 1u) << 31);
            float f; memcpy(&f, &b, 4); return f;
        }
        case 1:                             /* signed zero */
            return (bits & 1u) ? -0.0f : 0.0f;
        case 2:                             /* +/- infinity (or finite remap) */
            if (finite_only) return ldexpf(1.0f, (int)(bits % 200) - 100);
            return (bits & 1u) ? -INFINITY : INFINITY;
        case 3: {                           /* NaN (or finite remap) */
            if (finite_only) return ldexpf((float)(1 + (rng_next() % 999)), (int)(bits % 240) - 120);
            uint32_t b = 0x7FC00000u | (rng_next() & 0x003FFFFFu);
            float f; memcpy(&f, &b, 4); return f;
        }
        case 4: case 5: case 6: case 7: {   /* "geometry-like" scaled decimals; ldexpf, not powf */
            int e = (int)(bits % 9) - 4;
            return ldexpf((float)((int)(rng_next() % 2001) - 1000), e);
        }
        default: {                          /* general finite random float */
            uint32_t b = (rng_next() & 0x7F7FFFFFu) | ((bits & 1u) << 31);
            float f; memcpy(&f, &b, 4); return f;
        }
    }
}

static void fill_mesh(Mesh* m, int degenerate) {
    memset(m, 0, sizeof(*m));
    m->nv = 1 + (int)(rng_next() % MAXV);
    for (int i = 0; i < m->nv; i++) {
        m->v[i].x = stress_float();
        m->v[i].y = stress_float();
        m->v[i].z = stress_float();
    }
    int nf = (int)(rng_next() % (MAXF + 1));
    for (int i = 0; i < nf; i++) {
        int a = (int)(rng_next() % (unsigned)m->nv);
        int b = (int)(rng_next() % (unsigned)m->nv);
        int c = (int)(rng_next() % (unsigned)m->nv);
        if (!degenerate && m->nv >= 3) {   /* distinct indices; nv<3 can only be degenerate anyway */
            while (b == a) b = (b + 1) % m->nv;
            while (c == a || c == b) c = (c + 1) % m->nv;
        }
        m->f[i].a = a; m->f[i].b = b; m->f[i].c = c;
    }
    m->nf = nf;
}

static uint32_t bits_of(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }

/* returns 0 on identical (per policy above), otherwise a bitmask of bad axes */
static int component_compare(float a, float b) {
    if (isnan(a) || isnan(b))
        return (isnan(a) && isnan(b)) ? 0 : 1;      /* NaN must stay NaN, never become finite */
    return bits_of(a) == bits_of(b) ? 0 : 1;
}

static int outputs_equal(const V3* a, const V3* b, int nv) {
    for (int i = 0; i < nv; i++) {
        if (component_compare(a[i].x, b[i].x) ||
            component_compare(a[i].y, b[i].y) ||
            component_compare(a[i].z, b[i].z))
            return 0;
    }
    return 1;
}

static int run_phase(int ntrials, int* finite_trials, int* nan_trials) {
    int failures = 0;
    for (int t = 0; t < ntrials; t++) {
        Mesh m;
        fill_mesh(&m, (int)(rng_next() & 1u));

        int input_has_nan = 0;
        for (int i = 0; i < m.nv; i++)
            if (isnan(m.v[i].x) || isnan(m.v[i].y) || isnan(m.v[i].z))
                input_has_nan = 1;

        /* poison both output buffers with the SAME random bit pattern so the
         * zero-initialization contract is compared too (only out[0..nv-1]
         * is defined to be touched). */
        V3 out_new[MAXV], out_ref[MAXV];
        for (int k = 0; k < MAXV; k++) {
            out_new[k] = out_ref[k] = (V3){ stress_float(), stress_float(), stress_float() };
        }

        int r_new = vertex_normals(&m, out_new);
        int r_ref = vertex_normals_ref(&m, out_ref);

        int ok = (r_new == r_ref) && outputs_equal(out_new, out_ref, m.nv);
        if (input_has_nan) (*nan_trials)++; else (*finite_trials)++;
        if (!ok) {
            printf("MISMATCH trial=%d nf=%d nv=%d ret new=%d ref=%d input_nan=%d finite_only=%d\n",
                   t, m.nf, m.nv, r_new, r_ref, input_has_nan, finite_only);
            if (++failures >= 10) break;
        }
    }
    return failures;
}

int main(void) {
    enum { NTRIALS = 200000 };
    int failures = 0;
    int finite_trials = 0, nan_trials = 0;

    /* phase 1: full float space including NaN/Inf/subnormal stress */
    failures += run_phase(NTRIALS, &finite_trials, &nan_trials);

    /* phase 2: guaranteed-finite inputs, so EVERY comparison is exact-bit */
    int f2 = 0, n2 = 0;
    finite_only = 1;
    failures += run_phase(NTRIALS, &f2, &n2);
    finite_only = 0;

    /* targeted: tiny coordinates whose accumulated length lands below 1e-12 */
    for (int t = 0; t < 1000 && failures < 10; t++) {
        Mesh m; memset(&m, 0, sizeof(m));
        m.nv = 3;
        float s = ldexpf(1.0f, -70 - (int)(rng_next() % 40));   /* 2^-70 .. 2^-109 scale */
        for (int i = 0; i < 3; i++) { m.v[i].x = stress_float() * s;
                                      m.v[i].y = stress_float() * s;
                                      m.v[i].z = stress_float() * s; }
        m.f[0] = (Tri){0,1,2}; m.nf = 1;
        V3 a[MAXV], b[MAXV];
        memset(a, 0x5A, sizeof(a)); memset(b, 0x5A, sizeof(b));
        int ra = vertex_normals(&m, a), rb = vertex_normals_ref(&m, b);
        if (ra != rb || !outputs_equal(a, b, 3)) {
            printf("MISMATCH tiny trial=%d ret %d/%d\n", t, ra, rb);
            if (++failures >= 10) break;
        }
    }

    if (failures == 0)
        printf("bit_identity: %d mixed trials (%d finite exact-bit, %d NaN-input NaN-ness) + "
               "%d guaranteed-finite exact-bit trials + 1000 tiny trials, all equivalent\n",
               NTRIALS, finite_trials, nan_trials, f2);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
