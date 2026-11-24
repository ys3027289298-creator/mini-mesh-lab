
/* unit_tests: branch-level tests for the four split units U1..U4.
 * These exercise branches the original end-to-end main() cannot reach. */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../src/normal.h"

V3 face_normal(const Mesh* m, int i);
void accumulate_face(const Mesh* m, int i, V3 n, V3* out);
int normalize_vertex(V3* v);

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static uint32_t bits_of(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }
static int same_bits(float a, float b) { return bits_of(a) == bits_of(b); }

int main(void) {
    /* ================= U1: face_normal ================= */
    Mesh m; memset(&m, 0, sizeof(m));
    m.nv = 6;
    m.v[0] = (V3){0,0,0}; m.v[1] = (V3){1,0,0}; m.v[2] = (V3){0,1,0};
    m.v[3] = (V3){0,0,1};
    m.v[4] = (V3){5,5,5};                       /* translated copy of tri 0 */
    m.v[5] = (V3){5,5,5};                       /* repeated vertex -> degenerate */
    m.f[0] = (Tri){0,1,2};                      /* CCW: normal +z, (0,0,1) */
    m.f[1] = (Tri){0,2,1};                      /* reversed winding: (0,0,-1) */
    m.f[2] = (Tri){0,1,3};                      /* cross((1,0,0),(0,0,1)) = (0,-1,0) */
    m.f[3] = (Tri){4,5,5};                      /* zero-area: cross must be (0,0,0) exactly */
    m.f[4] = (Tri){4,4,4};                      /* all-same-vertex */
    m.nf = 5;

    V3 n0 = face_normal(&m, 0);
    check("u1_basic_plus_z", same_bits(n0.x,0.0f) && same_bits(n0.y,0.0f) && same_bits(n0.z,1.0f));
    V3 n1 = face_normal(&m, 1);
    check("u1_winding_sign", same_bits(n1.z,-1.0f));
    V3 n2 = face_normal(&m, 2);
    check("u1_second_axis", same_bits(n2.x,0.0f) && same_bits(n2.y,-1.0f) && same_bits(n2.z,0.0f));
    V3 n3 = face_normal(&m, 3);
    check("u1_degenerate_zero", same_bits(n3.x,0.0f) && same_bits(n3.y,0.0f) && same_bits(n3.z,0.0f));
    V3 n4 = face_normal(&m, 4);
    check("u1_all_same_vertex", same_bits(n4.x,0.0f) && same_bits(n4.y,0.0f) && same_bits(n4.z,0.0f));
    /* translation invariance: tri shifted by (5,5,5) must give identical bits */
    Mesh mt; memset(&mt, 0, sizeof(mt));
    mt.nv = 3;
    mt.v[0] = (V3){5,5,5}; mt.v[1] = (V3){6,5,5}; mt.v[2] = (V3){5,6,5};
    mt.f[0] = (Tri){0,1,2};
    V3 nt = face_normal(&mt, 0);
    check("u1_translation_invariant", same_bits(nt.x,n0.x) && same_bits(nt.y,n0.y) && same_bits(nt.z,n0.z));

    /* ================= U2: accumulate_face ================= */
    V3 acc[MAXV];
    for (int i = 0; i < MAXV; i++) acc[i] = (V3){0,0,0};
    V3 fn = face_normal(&m, 0);               /* (0,0,1) */
    accumulate_face(&m, 0, fn, acc);
    check("u2_scatter_corners",
          same_bits(acc[0].z,1.0f) && same_bits(acc[1].z,1.0f) && same_bits(acc[2].z,1.0f) &&
          same_bits(acc[0].x,0.0f) && same_bits(acc[0].y,0.0f));
    check("u2_untouched_vertex",
          same_bits(acc[3].x,0.0f) && same_bits(acc[3].y,0.0f) && same_bits(acc[3].z,0.0f));
    /* repeated index: one corner accumulates twice — an ordering/aliasing case
     * the e2e meshes never build. */
    V3 acc2[MAXV]; memset(acc2, 0, sizeof(acc2));
    Tri saved = m.f[0];
    m.f[0] = (Tri){0,0,1};
    V3 fn2 = face_normal(&m, 0);              /* b-a == 0 -> zero vector */
    accumulate_face(&m, 0, fn2, acc2);
    check("u2_repeated_index_double_add",
          same_bits(acc2[0].x,0.0f) && same_bits(acc2[0].y,0.0f) && same_bits(acc2[0].z,0.0f) &&
          same_bits(acc2[1].x,0.0f) && same_bits(acc2[1].y,0.0f) && same_bits(acc2[1].z,0.0f));
    /* in-place accumulation order: two sequential adds must equal the literal
     * sequential float sum (protects against a temp-buffer "refactor"). */
    V3 acc3[MAXV]; memset(acc3, 0, sizeof(acc3));
    m.f[0] = saved;
    V3 q = (V3){1.0f, 1e7f, -3.25f};
    accumulate_face(&m, 0, q, acc3);
    accumulate_face(&m, 0, q, acc3);
    float expect_y = ((0.0f + q.y) + q.y);    /* exact evaluation order */
    check("u2_inplace_add_order", same_bits(acc3[0].y, expect_y) && same_bits(acc3[0].x, 2.0f));

    /* ================= U3: normalize_vertex ================= */
    V3 u = (V3){3.0f, 0.0f, 0.0f};
    check("u3_normalize_ok", normalize_vertex(&u) == 1 &&
          same_bits(u.x, 1.0f) && same_bits(u.y, 0.0f) && same_bits(u.z, 0.0f));

    V3 z = (V3){0.0f, 0.0f, 0.0f};            /* positive zero */
    check("u3_zero_rejected", normalize_vertex(&z) == 0 &&
          same_bits(z.x,0.0f) && same_bits(z.y,0.0f) && same_bits(z.z,0.0f));

    V3 nz; memset(&nz, 0, sizeof(nz));
    nz.x = -0.0f;                              /* signed zero must be preserved, not overwritten */
    int rnz = normalize_vertex(&nz);
    check("u3_signed_zero_preserved", rnz == 0 && (bits_of(nz.x) == bits_of(-0.0f)));

    /* sub-threshold NONZERO vector: l in (0, 1e-12). Must return 0 and leave
     * every component bit-untouched. */
    V3 tiny = (V3){1e-20f, 0.0f, 0.0f};
    uint32_t tx = bits_of(tiny.x);
    check("u3_subthreshold_untouched", normalize_vertex(&tiny) == 0 && bits_of(tiny.x) == tx);

    /* just above threshold: normalized to unit */
    V3 near_hi = (V3){2e-12f, 0.0f, 0.0f};
    check("u3_above_threshold", normalize_vertex(&near_hi) == 1 &&
          fabsf(near_hi.x - 1.0f) < 1e-6f);

    /* NaN input: len is NaN, every comparison is false -> treated as
     * normalizable and divided; the result stays NaN. Payload bits are not
     * guaranteed by the C standard (hardware may rewrite NaN/NaN), so only
     * the branch outcome and NaN-ness are asserted; exact bit behavior vs.
     * the original monolith is covered by bit_identity. */
    V3 nv; memcpy(&nv.x, &(uint32_t){0x7FC01234u}, 4); nv.y = 0.0f; nv.z = 0.0f;
    int rnan = normalize_vertex(&nv);
    check("u3_nan_branch", rnan == 1 && isnan(nv.x));

    /* ================= U4: vertex_normals ================= */
    Mesh empty; memset(&empty, 0, sizeof(empty));
    empty.nv = 1;
    V3 dummy[MAXV];
    check("u4_no_faces", vertex_normals(&empty, dummy) == -1);

    /* all-degenerate mesh: every accumulator stays zero -> bad == nv */
    Mesh deg; memset(&deg, 0, sizeof(deg));
    deg.nv = 3;
    deg.v[0] = (V3){1,1,1}; deg.v[1] = (V3){1,1,1}; deg.v[2] = (V3){1,1,1};
    deg.f[0] = (Tri){0,1,2}; deg.nf = 1;
    V3 dn[MAXV];
    int bad = vertex_normals(&deg, dn);
    check("u4_all_bad_count", bad == 3);
    check("u4_bad_left_zero", same_bits(dn[0].x,0.0f) && same_bits(dn[1].z,0.0f));

    /* partial: one well-formed tri plus a degenerate one sharing no vertex */
    Mesh mix; memset(&mix, 0, sizeof(mix));
    mix.nv = 6;
    mix.v[0] = (V3){0,0,0}; mix.v[1] = (V3){1,0,0}; mix.v[2] = (V3){0,1,0};
    mix.v[3] = mix.v[4] = mix.v[5] = (V3){9,9,9};
    mix.f[0] = (Tri){0,1,2};
    mix.f[1] = (Tri){3,4,5};
    mix.nf = 2;
    V3 mn[MAXV];
    int mbad = vertex_normals(&mix, mn);
    check("u4_partial_bad", mbad == 3);
    check("u4_good_still_unit", fabsf(sqrtf(mn[0].x*mn[0].x + mn[0].y*mn[0].y + mn[0].z*mn[0].z) - 1.0f) < 1e-6f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
