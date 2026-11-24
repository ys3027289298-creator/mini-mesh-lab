
/* lightmap: UV-island padding for lightmap baking.
 * Padding width is derived from texture size (filter footprint) and
 * island spacing (gap between islands), not a fixed texel count.
 *
 * Model: two rectangular islands on one N x N lightmap.
 *   island A bakes to 1.0 (lit), island B to 0.0 (dark), background 0.0.
 * Bleed: after a box filter of radius FILT_R (simulates bilinear/mip
 * sampling), island-A interior texels pulled below THRESH are counted.
 *
 * coupled  mode: p = min(FILT_R, gap/2), multi-source BFS fills background
 *                only, so each island's ring is capped by the actual gap.
 * fixed    mode: p = FILT_R always, stamped from the original silhouette
 *                with no clipping -- the texel-decoupled bug, for contrast.
 */
#include <stdio.h>
#include <string.h>

#define N        128
#define FILT_R   (N / 32)   /* filter radius in texels; scales with texture size */
#define THRESH   0.95f

#define IW 40
#define IH 64
#define AX0 4
#define AY0 32

enum { SRC_BG = -1, SRC_A = 0, SRC_B = 1, SRC_PADA = 2, SRC_PADB = 3 };

typedef struct {
    float c[N * N];
    int   src[N * N];
    unsigned char interiorA[N * N]; /* island A interior mask, fixed at rasterize */
} LM;

static void rasterize(LM* lm, int gap) {
    for (int i = 0; i < N * N; i++) { lm->c[i] = 0.0f; lm->src[i] = SRC_BG; lm->interiorA[i] = 0; }
    for (int y = AY0; y < AY0 + IH; y++)
        for (int x = AX0; x < AX0 + IW; x++) {
            lm->c[y * N + x] = 1.0f; lm->src[y * N + x] = SRC_A; lm->interiorA[y * N + x] = 1;
        }
    int bx0 = AX0 + IW + gap;
    for (int y = AY0; y < AY0 + IH; y++)
        for (int x = bx0; x < bx0 + IW; x++) {
            lm->c[y * N + x] = 0.0f; lm->src[y * N + x] = SRC_B;
        }
}

/* coupled: multi-source Chebyshev BFS from both islands, background only,
 * depth <= p. Where the gap is narrower than 2p the two rings meet in the
 * middle instead of overlapping, so the effective width is min(p, gap/2). */
static void dilate_coupled(LM* lm, int p) {
    static int q[N * N], dist[N * N];
    static float col[N * N];
    int qh = 0, qt = 0;
    for (int i = 0; i < N * N; i++)
        if (lm->src[i] == SRC_A || lm->src[i] == SRC_B) {
            dist[i] = 0; col[i] = lm->c[i]; q[qt++] = i;
        }
    while (qh < qt) {
        int i = q[qh++];
        if (dist[i] >= p) continue;
        int x = i % N, y = i / N;
        int isA = (lm->src[i] == SRC_A || lm->src[i] == SRC_PADA);
        for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
                if (!dx && !dy) continue;
                int nx = x + dx, ny = y + dy;
                if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
                int j = ny * N + nx;
                if (lm->src[j] != SRC_BG) continue;
                lm->src[j] = isA ? SRC_PADA : SRC_PADB;
                lm->c[j] = col[i];
                dist[j] = dist[i] + 1;
                col[j] = col[i];
                q[qt++] = j;
            }
    }
}

/* decoupled (buggy): fixed p texels per island, distance measured from the
 * original interior snapshot, stamped over whatever is already there --
 * including the neighbouring island's interior. */
static void dilate_fixed(LM* lm, int p) {
    static unsigned char interior[2][N * N];
    static int dist[N * N], q[N * N];
    static unsigned char vis[N * N];
    for (int i = 0; i < N * N; i++) {
        interior[0][i] = (lm->src[i] == SRC_A);
        interior[1][i] = (lm->src[i] == SRC_B);
    }
    for (int island = 0; island < 2; island++) {
        int pad = island ? SRC_PADB : SRC_PADA;
        float col = island ? 0.0f : 1.0f;
        memset(vis, 0, sizeof(vis));
        int qh = 0, qt = 0;
        for (int i = 0; i < N * N; i++)
            if (interior[island][i]) { vis[i] = 1; dist[i] = 0; q[qt++] = i; }
        while (qh < qt) {
            int i = q[qh++];
            if (dist[i] >= p) continue;
            int x = i % N, y = i / N;
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    if (!dx && !dy) continue;
                    int nx = x + dx, ny = y + dy;
                    if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
                    int j = ny * N + nx;
                    if (vis[j]) continue;
                    vis[j] = 1; dist[j] = dist[i] + 1; q[qt++] = j;
                }
        }
        for (int k = 0; k < qt; k++) {
            int i = q[k];
            if (interior[island][i]) continue;
            lm->src[i] = pad;
            lm->c[i] = col;
        }
    }
}

static void box_filter(const float* in, float* out) {
    static float tmp[N * N];
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            float s = 0.0f; int n = 0;
            for (int dx = -FILT_R; dx <= FILT_R; dx++) {
                int xx = x + dx;
                if (xx < 0 || xx >= N) continue;
                s += in[y * N + xx]; n++;
            }
            tmp[y * N + x] = s / n;
        }
    for (int y = 0; y < N; y++)
        for (int x = 0; x < N; x++) {
            float s = 0.0f; int n = 0;
            for (int dy = -FILT_R; dy <= FILT_R; dy++) {
                int yy = y + dy;
                if (yy < 0 || yy >= N) continue;
                s += tmp[yy * N + x]; n++;
            }
            out[y * N + x] = s / n;
        }
}

/* island-A interior texels pulled below THRESH after filtering */
static int bleed_count(const LM* lm) {
    static float out[N * N];
    box_filter(lm->c, out);
    int n = 0;
    for (int i = 0; i < N * N; i++)
        if (lm->interiorA[i] && out[i] < THRESH) n++;
    return n;
}

/* min over island-A interior texels of (Chebyshev distance to the nearest
 * texel not coloured by A) - 1 == width of the safe ring at its narrowest */
static int min_margin(const LM* lm) {
    static int dist[N * N], q[N * N];
    int qh = 0, qt = 0;
    for (int i = 0; i < N * N; i++) {
        if (lm->src[i] != SRC_A && lm->src[i] != SRC_PADA) { dist[i] = 0; q[qt++] = i; }
        else dist[i] = -1;
    }
    while (qh < qt) {
        int i = q[qh++];
        int x = i % N, y = i / N;
        for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
                if (!dx && !dy) continue;
                int nx = x + dx, ny = y + dy;
                if (nx < 0 || ny < 0 || nx >= N || ny >= N) continue;
                int j = ny * N + nx;
                if (dist[j] < 0) { dist[j] = dist[i] + 1; q[qt++] = j; }
            }
    }
    int m = 1 << 30;
    for (int i = 0; i < N * N; i++)
        if (lm->interiorA[i] && dist[i] >= 0 && dist[i] - 1 < m) m = dist[i] - 1;
    return m;
}

static double pad_occupancy(const LM* lm) {
    int n = 0;
    for (int i = 0; i < N * N; i++)
        if (lm->src[i] == SRC_PADA || lm->src[i] == SRC_PADB) n++;
    return 100.0 * n / (N * N);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    static LM lm;
    const int gaps[3] = {2, 8, 32};
    int bleed_base[3], bleed_c[3], bleed_f[3], margin_c[3], margin_f[3];
    double occ_c[3];

    printf("# coupled padding: p = min(filt_r=%d, gap/2) texels, N=%d\n", FILT_R, N);
    for (int k = 0; k < 3; k++) {
        int gap = gaps[k];
        rasterize(&lm, gap);
        bleed_base[k] = bleed_count(&lm);
        rasterize(&lm, gap);
        dilate_coupled(&lm, FILT_R);
        bleed_c[k] = bleed_count(&lm);
        occ_c[k] = pad_occupancy(&lm);
        margin_c[k] = min_margin(&lm);
        printf("coupled gap=%2d bleed=%4d pad_occ=%5.2f%% margin=%2d tex (%.4f uv) [baseline bleed=%d]\n",
               gap, bleed_c[k], occ_c[k], margin_c[k], margin_c[k] / (double)N, bleed_base[k]);
    }

    printf("# decoupled fixed-texel padding: p = %d tex (%.4f uv) at every gap\n",
           FILT_R, FILT_R / (double)N);
    for (int k = 0; k < 3; k++) {
        rasterize(&lm, gaps[k]);
        dilate_fixed(&lm, FILT_R);
        bleed_f[k] = bleed_count(&lm);
        margin_f[k] = min_margin(&lm);
        printf("fixed   gap=%2d bleed=%4d pad_w=%d tex (%.4f uv) margin=%d tex\n",
               gaps[k], bleed_f[k], FILT_R, FILT_R / (double)N, margin_f[k]);
    }
    printf("# gap/2 allows pad <= %d/%d/%d tex = %.4f/%.4f/%.4f uv\n",
           gaps[0] / 2, gaps[1] / 2, gaps[2] / 2,
           (gaps[0] / 2) / (double)N, (gaps[1] / 2) / (double)N, (gaps[2] / 2) / (double)N);

    check("narrow_gap_bleed_reduced", bleed_c[0] < bleed_base[0]);
    check("wide_gaps_clean", bleed_c[0+1] == 0 && bleed_c[2] == 0);
    check("margin_tracks_gap", margin_c[0] == 1 && margin_c[1] == FILT_R && margin_c[2] == FILT_R);
    check("fixed_mode_fails_narrow", bleed_f[0] > bleed_c[0]);
    check("fixed_mode_overlaps_island", margin_f[0] < 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
