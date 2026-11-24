
/* bench_ext: beyond-2000 points to validate the frame-budget extrapolation. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "collision.h"

#define BLOCKS 5
#define WARMUP 1
#define STEP_FRAC 0.01f

static unsigned g_seed = 0x9e3779b9u;
static float rnd(void) {
    unsigned x = g_seed;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    g_seed = x;
    return (x >> 8) * (1.0f / 16777216.0f);
}
static long long now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long long)ts.tv_sec * 1000000000LL + ts.tv_nsec;
}
static int cmp_ll(const void* a, const void* b) {
    long long x = *(const long long*)a, y = *(const long long*)b;
    return (x > y) - (x < y);
}

int main(void) {
    static const int ns[] = { 2500, 3000, 3500, 4000 };
    int nmax = 4000;
    Body* bodies = malloc(sizeof(Body) * nmax);
    Proxy* px = malloc(sizeof(Proxy) * nmax);
    Pair* pairs = malloc(sizeof(Pair) * (size_t)nmax * (nmax - 1) / 2);
    float* verts = malloc(sizeof(float) * nmax * BODY_VERTS * 3);
    if (!bodies || !px || !pairs || !verts) { fprintf(stderr, "oom\n"); return 1; }

    printf("n,frames,total_ns,proxy_ns,broad_ns,pairs\n");
    for (size_t t = 0; t < sizeof(ns)/sizeof(ns[0]); t++) {
        int n = ns[t];
        float L = cbrtf((float)n), step = L * STEP_FRAC;
        int cap = (int)((size_t)n * (n - 1) / 2);
        for (int i = 0; i < n; i++) {
            bodies[i].base = i * BODY_VERTS;
            body_make_mesh(verts, &bodies[i], &g_seed);
            bodies[i].q[0] = 1.0f; bodies[i].q[1] = bodies[i].q[2] = bodies[i].q[3] = 0.0f;
            for (int k = 0; k < 3; k++) bodies[i].pos[k] = rnd() * L;
        }
        int frames = 3;
        long long btot[BLOCKS], bprx[BLOCKS], bbrd[BLOCKS];
        int last_pairs = 0;
        volatile unsigned sink = 0;
        for (int blk = 0; blk < WARMUP + BLOCKS; blk++) {
            int timed = blk >= WARMUP, bi = blk - WARMUP;
            long long t0, t1; unsigned cs = 0;
            t0 = now_ns();
            for (int f = 0; f < frames; f++) proxy_build(verts, bodies, px, n);
            t1 = now_ns(); if (timed) bprx[bi] = (t1-t0)/frames;
            t0 = now_ns();
            for (int f = 0; f < frames; f++) last_pairs = broadphase(px, n, pairs, cap);
            t1 = now_ns(); if (timed) bbrd[bi] = (t1-t0)/frames;
            t0 = now_ns();
            for (int f = 0; f < frames; f++) {
                scene_advance(bodies, n, L, step, &g_seed);
                proxy_build(verts, bodies, px, n);
                int np = broadphase(px, n, pairs, cap);
                cs += np;
            }
            t1 = now_ns(); if (timed) btot[bi] = (t1-t0)/frames;
            sink ^= cs;
        }
        qsort(btot, BLOCKS, sizeof(long long), cmp_ll);
        qsort(bprx, BLOCKS, sizeof(long long), cmp_ll);
        qsort(bbrd, BLOCKS, sizeof(long long), cmp_ll);
        printf("%d,%d,%lld,%lld,%lld,%d\n", n, frames,
               btot[BLOCKS/2], bprx[BLOCKS/2], bbrd[BLOCKS/2], last_pairs);
        fflush(stdout);
        fprintf(stderr, "n=%d total=%.4f ms proxy=%.4f broad=%.4f (checksum %x)\n",
                n, btot[BLOCKS/2]/1e6, bprx[BLOCKS/2]/1e6, bbrd[BLOCKS/2]/1e6, sink);
    }
    free(bodies); free(px); free(pairs); free(verts);
    return 0;
}
