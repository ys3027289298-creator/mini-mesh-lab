
/* collision: per-frame proxy generation, brute broadphase, sphere narrowphase. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "collision.h"

/* deterministic xorshift32, returns float in [0,1) */
static float rnd(unsigned* s) {
    unsigned x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return (x >> 8) * (1.0f / 16777216.0f);
}

#define MESH_SCALE 0.03f   /* local point cloud lives in [-MESH_SCALE, +MESH_SCALE]^3 */

void body_make_mesh(float* verts, Body* body, unsigned* rng) {
    float r2max = 0.0f;
    for (int v = 0; v < BODY_VERTS; v++) {
        float* p = &verts[(body->base + v) * 3];
        for (int k = 0; k < 3; k++) p[k] = (rnd(rng) * 2.0f - 1.0f) * MESH_SCALE;
        float r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        if (r2 > r2max) r2max = r2;
    }
    body->r = sqrtf(r2max);
}

/* quaternion rotate vector p by unit quaternion q (qw,qx,qy,qz) */
static void qrot(const float q[4], const float p[3], float out[3]) {
    float w = q[0], x = q[1], y = q[2], z = q[3];
    /* v' = p + 2 * cross(q.xyz, cross(q.xyz, p) + w*p) */
    float u[3] = { x, y, z };
    float t[3];
    t[0] = 2.0f * (u[1] * p[2] - u[2] * p[1] + w * p[0]);
    t[1] = 2.0f * (u[2] * p[0] - u[0] * p[2] + w * p[1]);
    t[2] = 2.0f * (u[0] * p[1] - u[1] * p[0] + w * p[2]);
    out[0] = p[0] + (u[1] * t[2] - u[2] * t[1]);
    out[1] = p[1] + (u[2] * t[0] - u[0] * t[2]);
    out[2] = p[2] + (u[0] * t[1] - u[1] * t[0]);
}

void proxy_build(const float* verts, const Body* bodies, Proxy* proxies, int n) {
    for (int i = 0; i < n; i++) {
        const Body* b = &bodies[i];
        AABB* box = &proxies[i].box;
        box->lo[0] = box->lo[1] = box->lo[2] = 1e30f;
        box->hi[0] = box->hi[1] = box->hi[2] = -1e30f;
        for (int v = 0; v < BODY_VERTS; v++) {
            const float* p = &verts[(b->base + v) * 3];
            float w[3];
            qrot(b->q, p, w);
            if (w[0] < box->lo[0]) box->lo[0] = w[0];
            if (w[1] < box->lo[1]) box->lo[1] = w[1];
            if (w[2] < box->lo[2]) box->lo[2] = w[2];
            if (w[0] > box->hi[0]) box->hi[0] = w[0];
            if (w[1] > box->hi[1]) box->hi[1] = w[1];
            if (w[2] > box->hi[2]) box->hi[2] = w[2];
        }
        for (int k = 0; k < 3; k++) {
            box->lo[k] += b->pos[k];
            box->hi[k] += b->pos[k];
        }
    }
}

static int aabb_overlap(const AABB* x, const AABB* y) {
    for (int k = 0; k < 3; k++)
        if (x->hi[k] < y->lo[k] || y->hi[k] < x->lo[k]) return 0;
    return 1;
}

int broadphase(const Proxy* proxies, int n, Pair* out, int cap) {
    int cnt = 0;
    for (int i = 0; i < n;  i++) {
        for (int j = i + 1; j < n; j++) {
            if (!aabb_overlap(&proxies[i].box, &proxies[j].box)) continue;
            if (cnt < cap) out[cnt].a = i, out[cnt].b = j;
            cnt++;
        }
    }
    return cnt;
}

int narrowphase(const Body* bodies, const Pair* pairs, int npair) {
    int contacts = 0;
    for (int p = 0; p < npair; p++) {
        const Body* a = &bodies[pairs[p].a];
        const Body* b = &bodies[pairs[p].b];
        float rr = a->r + b->r;
        float d2 = 0.0f;
        for (int k = 0; k < 3; k++) {
            float d = a->pos[k] - b->pos[k];
            d2 += d * d;
        }
        if (d2 <= rr * rr) contacts++;
    }
    return contacts;
}

static void qnormalize(float q[4]) {
    float l = sqrtf(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    if (l > 0.0f) for (int k = 0; k < 4; k++) q[k] /= l;
}

static void qmul(const float a[4], const float b[4], float o[4]) {
    o[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    o[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
    o[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
    o[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
}

void scene_advance(Body* bodies, int n, float L, float step, unsigned* rng) {
    for (int i = 0; i < n; i++) {
        Body* b = &bodies[i];
        for (int k = 0; k < 3; k++) {
            b->pos[k] += (rnd(rng) * 2.0f - 1.0f) * step;
            while (b->pos[k] < 0.0f)  b->pos[k] += L;
            while (b->pos[k] >= L)    b->pos[k] -= L;
        }
        float axis[3] = { rnd(rng) * 2.0f - 1.0f,
                          rnd(rng) * 2.0f - 1.0f,
                          rnd(rng) * 2.0f - 1.0f };
        float al = sqrtf(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
        float ang = rnd(rng) * 0.1f;
        float dq[4], out[4];
        if (al > 0.0f) {
            float s = sinf(ang * 0.5f) / al;
            dq[0] = cosf(ang * 0.5f); dq[1] = axis[0]*s;
            dq[2] = axis[1]*s;        dq[3] = axis[2]*s;
        } else {
            dq[0] = 1.0f; dq[1] = dq[2] = dq[3] = 0.0f;
        }
        qmul(dq, b->q, out);
        memcpy(b->q, out, sizeof(float) * 4);
        qnormalize(b->q);
    }
}

#ifndef COLLISION_NO_MAIN
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    enum { NB = 3, NV = NB * BODY_VERTS };
    static float verts[NV * 3];
    Body bodies[NB];
    unsigned seed = 7;
    memset(bodies, 0, sizeof(bodies));
    for (int i = 0; i < NB; i++) {
        bodies[i].base = i * BODY_VERTS;
        body_make_mesh(verts, &bodies[i], &seed);
        bodies[i].q[0] = 1.0f;   /* identity orientation */
        for (int k = 0; k < 3; k++) bodies[i].pos[k] = 0.0f;
    }

    /* two bodies near the origin, one far away. AABBs span about
       +/-MESH_SCALE, so a center gap of 0.05 keeps the boxes overlapping
       while the bounding spheres (r ~ 0.052) also overlap. */
    bodies[1].pos[0] = 0.05f;
    bodies[2].pos[0] = 100.0f;

    Proxy px[NB];
    proxy_build(verts, bodies, px, NB);

    /* identity rotation: proxy extents equal the raw mesh extents */
    float lo0 = 1e30f, hi0 = -1e30f;
    for (int v = 0; v < BODY_VERTS; v++) {
        float x = verts[v * 3];
        if (x < lo0) lo0 = x;
        if (x > hi0) hi0 = x;
    }
    check("proxy_identity_extent",
          fabsf(px[0].box.lo[0] - lo0) < 1e-5f &&
          fabsf(px[0].box.hi[0] - hi0) < 1e-5f);
    check("proxy_valid", px[0].box.lo[0] < px[0].box.hi[0]);

    /* rotation must not change the AABB size beyond numerical noise once a
       full 180deg turn is applied (symmetric box extents) */
    Body rotb = bodies[0];
    rotb.q[0] = 0.0f; rotb.q[1] = 1.0f;   /* 180 deg about x */
    Proxy rp;
    proxy_build(verts, &rotb, &rp, 1);
    float d_lo = fabsf(rp.box.lo[0] - px[0].box.lo[0]);
    check("proxy_rotation_x", d_lo < 1e-4f);

    Pair pairs[8];
    int np = broadphase(px, NB, pairs, 8);
    check("broadphase_one_pair", np == 1);
    check("broadphase_pair_idx", np == 1 && pairs[0].a == 0 && pairs[0].b == 1);

    /* cap=0 stores nothing but still reports the true count */
    check("cap_truncation", broadphase(px, NB, pairs, 0) == 1);

    int contacts = narrowphase(bodies, pairs, np);
    check("narrowphase_contact", contacts == 1);

    /* move body 1 far out of reach: no broadphase candidate, hence no contact */
    bodies[1].pos[0] = 1.0f;
    proxy_build(verts, bodies, px, NB);
    np = broadphase(px, NB, pairs, 8);
    check("separated_no_pair", np == 0);

    /* advance wraps every component back into [0,L) */
    scene_advance(bodies, NB, 5.0f, 0.05f, &seed);
    int wrapped = 1;
    for (int i = 0; i < NB; i++)
        for (int k = 0; k < 3; k++)
            if (bodies[i].pos[k] < 0.0f || bodies[i].pos[k] >= 5.0f) wrapped = 0;
    check("wrap_in_box", wrapped);

    /* quaternions stay unit length */
    int unit = 1;
    for (int i = 0; i < NB; i++) {
        float l2 = bodies[i].q[0]*bodies[i].q[0] + bodies[i].q[1]*bodies[i].q[1]
                 + bodies[i].q[2]*bodies[i].q[2] + bodies[i].q[3]*bodies[i].q[3];
        if (fabsf(l2 - 1.0f) > 1e-4f) unit = 0;
    }
    check("quaternion_unit", unit);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
