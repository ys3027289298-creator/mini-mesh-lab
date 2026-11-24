
/* collision: broadphase proxy generation + AABB broadphase + sphere narrowphase. */
#ifndef COLLISION_H
#define COLLISION_H

#define BODY_VERTS 256   /* local mesh vertices per object */

typedef struct { float lo[3], hi[3]; } AABB;
typedef struct { AABB box; } Proxy;               /* one broadphase proxy per object */
typedef struct { int a, b; } Pair;                /* candidate object pair */

typedef struct {
    float pos[3];      /* world translation */
    float q[4];        /* orientation quaternion: qw, qx, qy, qz */
    float r;           /* bounding radius of the local mesh */
    int   base;        /* first vertex in the shared pool */
} Body;

/* fill body's BODY_VERTS block in the pool with a random point cloud and
   set its bounding radius; body->base must already be assigned */
void body_make_mesh(float* verts, Body* body, unsigned* rng);

/* rebuild every proxy from its object mesh: O(n), BODY_VERTS transformed
   points per object (rotate + translate, then fit an AABB) */
void proxy_build(const float* verts, const Body* bodies, Proxy* proxies, int n);

/* brute-force AABB pair scan: O(n^2) overlap tests, appends candidate pairs.
   returns the true pair count; if cap is hit the output is truncated. */
int broadphase(const Proxy* proxies, int n, Pair* out, int cap);

/* exact sphere-sphere test over broadphase candidates: O(npair) */
int narrowphase(const Body* bodies, const Pair* pairs, int npair);

/* random-walk motion (translation + small rotation) with wrap into [0,L)^3:
   keeps scene density fixed */
void scene_advance(Body* bodies, int n, float L, float step, unsigned* rng);

#endif
