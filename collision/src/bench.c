
/* bench: stress the collision pipeline from 10 up to 2000 objects.
   each object carries a BODY_VERTS-vertex mesh; proxy generation transforms
   every vertex per frame. per phase we run BLOCKS blocks of `frames` frames
   each (frame count adapted so a block stays ~20-40ms); per-frame ns is the
   median over blocks. scene density is held fixed (L = cbrt(n)) so growth
   comes from algorithmic scaling, not crowding. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "collision.h"

#define BLOCKS      9
#define WARMUP      2
#define STEP_FRAC   0.01f
#define TARGET_BLOCK_NS 30000000LL   /* aim each timed block at ~30 ms */

static unsigned g_seed = 0x12345678u;
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

/* median per-frame ns over BLOCKS timed blocks */
static long long median_per_frame(const long long* block_ns, int frames) {
    long long pf[BLOCKS];
    for (int b = 0; b < BLOCKS; b++) pf[b] = block_ns[b] / frames;
    qsort(pf, BLOCKS, sizeof(long long), cmp_ll);
    return pf[BLOCKS / 2];
}

int main(void) {
    static const int ns[] = {
        10, 20, 40, 70, 100, 150, 220, 320, 460, 650, 900, 1250, 1700, 2000
    };
    int nmax = ns[sizeof(ns) / sizeof(ns[0]) - 1];

    Body*  bodies = malloc(sizeof(Body) * nmax);
    Proxy* px     = malloc(sizeof(Proxy) * nmax);
    Pair*  pairs  = malloc(sizeof(Pair) * (size_t)nmax * (nmax - 1) / 2);
    float* verts  = malloc(sizeof(float) * nmax * BODY_VERTS * 3);
    if (!bodies || !px || !pairs || !verts) { fprintf(stderr, "oom\n"); return 1; }

    volatile unsigned sink = 0;
    fprintf(stderr, "n, per-frame ms (median of %d blocks, mesh verts/object=%d)\n",
            BLOCKS, BODY_VERTS);
    printf("n,frames,total_ns,proxy_ns,broad_ns,narrow_ns,motion_ns,tests,pairs,contacts\n");

    for (size_t t = 0; t < sizeof(ns) / sizeof(ns[0]); t++) {
        int n = ns[t];
        float L = cbrtf((float)n);
        float step = L * STEP_FRAC;
        long long tests = (long long)n * (n - 1) / 2;
        int cap = (int)((size_t)n * (n - 1) / 2);

        for (int i = 0; i < n; i++) {
            bodies[i].base = i * BODY_VERTS;
            body_make_mesh(verts, &bodies[i], &g_seed);   /* sets bounding radius */
            bodies[i].q[0] = 1.0f; bodies[i].q[1] = bodies[i].q[2] = bodies[i].q[3] = 0.0f;
            for (int k = 0; k < 3; k++) bodies[i].pos[k] = rnd() * L;
        }

        /* adaptive frames: one cheap proxy frame costs ~n*BODY_VERTS work;
           scale the frame count so the block is long enough to time cleanly */
        int frames = (int)(TARGET_BLOCK_NS / (100LL * n * BODY_VERTS));
        if (frames < 4) frames = 4;
        if (frames > 4096) frames = 4096;

        long long b_proxy[BLOCKS], b_broad[BLOCKS], b_narrow[BLOCKS];
        long long b_motion[BLOCKS], b_total[BLOCKS];
        int last_pairs = 0, last_contacts = 0;

        for (int blk = 0; blk < WARMUP + BLOCKS; blk++) {
            int timed = blk >= WARMUP;
            int bi = blk - WARMUP;
            long long t0, t1;
            unsigned checksum = 0;

            /* phase 1: motion */
            t0 = now_ns();
            for (int f = 0; f < frames; f++)
                scene_advance(bodies, n, L, step, &g_seed);
            t1 = now_ns();
            if (timed) b_motion[bi] = t1 - t0;
            checksum += (unsigned)(bodies[0].pos[0] * 1000.0f);

            /* phase 2: proxy generation (mesh -> AABB) */
            t0 = now_ns();
            for (int f = 0; f < frames; f++) proxy_build(verts, bodies, px, n);
            t1 = now_ns();
            if (timed) b_proxy[bi] = t1 - t0;
            checksum += (unsigned)(px[n - 1].box.hi[0] * 1000.0f);

            /* phase 3: broadphase */
            t0 = now_ns();
            for (int f = 0; f < frames; f++)
                last_pairs = broadphase(px, n, pairs, cap);
            t1 = now_ns();
            if (timed) b_broad[bi] = t1 - t0;
            checksum += (unsigned)last_pairs;

            /* phase 4: narrowphase */
            t0 = now_ns();
            for (int f = 0; f < frames; f++)
                last_contacts = narrowphase(bodies, pairs, last_pairs);
            t1 = now_ns();
            if (timed) b_narrow[bi] = t1 - t0;
            checksum += (unsigned)last_contacts;

            /* full frame: everything one real frame does, timed together */
            t0 = now_ns();
            for (int f = 0; f < frames; f++) {
                scene_advance(bodies, n, L, step, &g_seed);
                proxy_build(verts, bodies, px, n);
                int np = broadphase(px, n, pairs, cap);
                last_contacts = narrowphase(bodies, pairs, np);
                checksum += np;
            }
            t1 = now_ns();
            if (timed) b_total[bi] = t1 - t0;
            sink ^= checksum;
        }

        long long pf_total  = median_per_frame(b_total, frames);
        long long pf_proxy  = median_per_frame(b_proxy, frames);
        long long pf_broad  = median_per_frame(b_broad, frames);
        long long pf_narrow = median_per_frame(b_narrow, frames);
        long long pf_motion = median_per_frame(b_motion, frames);

        printf("%d,%d,%lld,%lld,%lld,%lld,%lld,%lld,%d,%d\n",
               n, frames, pf_total, pf_proxy, pf_broad, pf_narrow, pf_motion,
               tests, last_pairs, last_contacts);
        fflush(stdout);
        fprintf(stderr, "n=%5d frames=%5d total=%9.4f ms proxy=%8.4f broad=%9.4f narrow=%7.5f pairs=%d\n",
                n, frames, pf_total / 1e6, pf_proxy / 1e6, pf_broad / 1e6,
                pf_narrow / 1e6, last_pairs);
    }

    fprintf(stderr, "checksum=%08x\n", sink);
    free(bodies); free(px); free(pairs); free(verts);
    return 0;
}
