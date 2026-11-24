/* mesh.c -- attribute-bank mesh for the redundancy experiment. */
#include "mesh.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

const int CH_COMP[CH_N] = { 3, 3, 2, 3 };
const char* const CH_NAME[CH_N] = { "pos", "nrm", "uv", "col" };

static float h01(unsigned a, unsigned b) {
    unsigned x = a * 73856093u ^ b * 19349663u ^ (b << 7);
    x ^= x >> 13; x *= 0x5bd1e995u; x ^= x >> 15;
    return (float)(x & 0xffffffu) / 16777216.0f;
}

static void ramp(float t, float* r, float* g, float* b) {
    *r = 0.35f + 0.65f * (0.5f + 0.5f * sinf(6.2831f * (t + 0.00f)));
    *g = 0.35f + 0.65f * (0.5f + 0.5f * sinf(6.2831f * (t + 0.33f)));
    *b = 0.35f + 0.65f * (0.5f + 0.5f * sinf(6.2831f * (t + 0.66f)));
}

/* surface sample at grid point (i,j) */
static void sample_at(int i, int j, int w, int h, V3* p, V3* n,
                      float* u, float* v, float* cr, float* cg, float* cb,
                      unsigned seed) {
    float a = (float)i / (float)(w - 1);
    float b = (float)j / (float)(h - 1);
    float x = (a - 0.5f) * 2.0f;
    float y = (b - 0.5f) * 2.0f;
    float z = 0.35f * sinf(2.0f * x) + 0.25f * cosf(1.5f * y);
    p->x = x; p->y = y; p->z = z;
    float dzdx = 0.70f * cosf(2.0f * x);
    float dzdy = -0.375f * sinf(1.5f * y);
    V3 t = { -dzdx, -dzdy, 1.0f };
    float l = sqrtf(t.x*t.x + t.y*t.y + t.z*t.z);
    n->x = t.x/l; n->y = t.y/l; n->z = t.z/l;
    *u = a; *v = b;
    float r, g, bl; ramp(a, &r, &g, &bl);
    float wob = 0.75f + 0.5f * h01((unsigned)((i * 131 + j * 37) * 7 + 3), seed);
    *cr = r * wob;
    *cg = g * (0.75f + 0.5f * h01((unsigned)((i * 131 + j * 37) * 7 + 4), seed));
    *cb = bl * (0.75f + 0.5f * h01((unsigned)((i * 131 + j * 37) * 7 + 5), seed));
}

static void shuffle(int* a, int n, unsigned seed) {
    for (int k = n - 1; k > 0; k--) {
        int r = (int)(h01((unsigned)(k + 1), seed) * (float)(k + 1));
        if (r > k) r = k;
        if (r < 0) r = 0;
        int t = a[k]; a[k] = a[r]; a[r] = t;
    }
}

static int grid_corner(int f, int k, int w) {
    int j = f / (2 * (w - 1));
    int r = f % (2 * (w - 1));
    int i = r / 2, s = r % 2;
    int a = j * w + i, b = j * w + i + 1;
    int c = (j + 1) * w + i + 1, d = (j + 1) * w + i;
    int tri[2][3] = { { a, b, c }, { a, c, d } };
    return tri[s][k];
}

static int rec_of_vertex(int v, int w, int h, int rx, int ry) {
    int i = v % w, j = v / w;
    int ri = (int)((long)i * rx / w);
    int rj = (int)((long)j * ry / h);
    if (ri >= rx) ri = rx - 1;
    if (rj >= ry) rj = ry - 1;
    return rj * rx + ri;
}

/* Build a w x h grid mesh.
 *
 * Position channel: bank = nv, stream references every vertex 1:1.
 * Attribute channel ch: bank = nv + pad[ch].  Its stream references
 * `live_n[ch]` slots, the deduplicated records.  The live slots are a
 * shuffled subset of [0, nv); every other slot -- the holes inside [0,nv)
 * plus the whole padded tail -- is stale data left by an earlier export and
 * is referenced by nothing.
 *
 * Because live slots are shuffled, attribute slot ids are not vertex ids.
 * Counting a channel through the position stream therefore reports all nv
 * low ids as used, which is wrong whenever live_n < nv. */
Mesh mesh_make(const char* name, int w, int h,
               const int div_x[CH_N], const int div_y[CH_N],
               const int bank_req[CH_N], float stale_eps, unsigned seed) {
    Mesh m; memset(&m, 0, sizeof(m));
    m.name = name;
    m.nv = w * h;
    m.nf = (w - 1) * (h - 1) * 2;
    int rx[CH_N], ry[CH_N];
    for (int ch = 0; ch < CH_N; ch++) {
        if (ch == CH_POS) {
            m.bank[ch] = m.nv;
            m.live_n[ch] = m.nv;
            rx[ch] = w; ry[ch] = h;
            continue;
        }
        int dx = div_x[ch] < 1 ? 1 : div_x[ch];
        int dy = div_y[ch] < 1 ? 1 : div_y[ch];
        rx[ch] = (w + dx - 1) / dx;
        ry[ch] = (h + dy - 1) / dy;
        if (rx[ch] < 1) rx[ch] = 1;
        if (ry[ch] < 1) ry[ch] = 1;
        m.live_n[ch] = rx[ch] * ry[ch];
        m.bank[ch] = bank_req[ch] < m.live_n[ch] ? m.live_n[ch] : bank_req[ch];
        if (m.bank[ch] > 65535) m.bank[ch] = 65535;
    }
    m.pos = calloc((size_t)m.bank[CH_POS], sizeof(V3));
    m.nrm = calloc((size_t)m.bank[CH_NRM], sizeof(V3));
    m.uv  = calloc((size_t)m.bank[CH_UV] * 2, sizeof(float));
    m.col = calloc((size_t)m.bank[CH_COL] * 3, sizeof(float));
    for (int ch = 0; ch < CH_N; ch++)
        m.idx[ch] = malloc((size_t)m.nf * 3 * sizeof(uint16_t));

    V3* full_p = malloc((size_t)m.nv * sizeof(V3));
    V3* full_n = malloc((size_t)m.nv * sizeof(V3));
    float* full_uv = malloc((size_t)m.nv * 2 * sizeof(float));
    float* full_c = malloc((size_t)m.nv * 3 * sizeof(float));
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            int v = j * w + i;
            V3 p, n; float u, t, r, g, b;
            sample_at(i, j, w, h, &p, &n, &u, &t, &r, &g, &b, seed);
            full_p[v] = p; full_n[v] = n;
            full_uv[2*v] = u; full_uv[2*v+1] = t;
            full_c[3*v+0] = r; full_c[3*v+1] = g; full_c[3*v+2] = b;
        }
    memcpy(m.pos, full_p, (size_t)m.nv * sizeof(V3));

    for (int ch = CH_NRM; ch <= CH_COL; ch++) {
        int nrec = m.live_n[ch];
        int bank = m.bank[ch];
        /* Live slots are drawn from [0, nv) only.  A row inside the vertex-id
         * range can therefore still be a dead hole, and the whole tail
         * [nv, bank) is dead.  A tool that walks vertex ids 0..nv-1 and
         * assumes "row i belongs to vertex i" sees every row it probes as
         * used, so it finds nothing to drop even when most of the bank is
         * dead.  That is the wrong-criterion failure this experiment shows. */
        int pull = m.nv < bank ? m.nv : bank;   /* live slots are < min(nv,bank) */
        int* all = malloc((size_t)pull * sizeof(int));
        for (int k = 0; k < pull; k++) all[k] = k;
        shuffle(all, pull, seed * 31u + (unsigned)ch * 1009u);
        int* slot_of_rec = all;
        int* rec_of_slot = malloc((size_t)m.bank[ch] * sizeof(int));
        for (int s = 0; s < m.bank[ch]; s++) rec_of_slot[s] = -1;
        for (int r = 0; r < nrec; r++) rec_of_slot[slot_of_rec[r]] = r;

        for (int r = 0; r < nrec; r++) {
            int ri = r % rx[ch], rj = r / rx[ch];
            int ci = (int)(((long)ri * 2 + 1) * w / (2 * rx[ch]));
            int cj = (int)(((long)rj * 2 + 1) * h / (2 * ry[ch]));
            if (ci >= w) ci = w - 1;
            if (cj >= h) cj = h - 1;
            int src = cj * w + ci;
            int sl = slot_of_rec[r];
            if (ch == CH_NRM) m.nrm[sl] = full_n[src];
            else if (ch == CH_UV) { m.uv[2*sl] = full_uv[2*src]; m.uv[2*sl+1] = full_uv[2*src+1]; }
            else {
                m.col[3*sl+0] = full_c[3*src+0];
                m.col[3*sl+1] = full_c[3*src+1];
                m.col[3*sl+2] = full_c[3*src+2];
            }
        }
        for (int sl = 0; sl < m.bank[ch]; sl++) {
            if (rec_of_slot[sl] >= 0) continue;
            int r = (int)(h01((unsigned)(sl * 13 + 5), seed * 7u + (unsigned)ch) *
                          (float)nrec);
            if (r >= nrec) r = nrec - 1;
            int rs = slot_of_rec[r];
            float e = stale_eps * (0.5f + h01((unsigned)(sl * 17 + 9), seed));
            if (ch == CH_NRM) {
                m.nrm[sl] = m.nrm[rs];
                m.nrm[sl].x += e * 0.01f; m.nrm[sl].z -= e * 0.02f;
            } else if (ch == CH_UV) {
                m.uv[2*sl] = m.uv[2*rs] + e * 0.01f;
                m.uv[2*sl+1] = m.uv[2*rs+1] - e * 0.01f;
            } else {
                m.col[3*sl+0] = m.col[3*rs+0] + e * 0.01f;
                m.col[3*sl+1] = m.col[3*rs+1];
                m.col[3*sl+2] = m.col[3*rs+2] - e * 0.01f;
            }
        }
        for (int f = 0; f < m.nf; f++)
            for (int k = 0; k < 3; k++) {
                int v = grid_corner(f, k, w);
                m.idx[ch][3*f+k] = (uint16_t)slot_of_rec[
                    rec_of_vertex(v, w, h, rx[ch], ry[ch])];
            }
        free(all);
        free(rec_of_slot);
    }
    for (int f = 0; f < m.nf; f++)
        for (int k = 0; k < 3; k++)
            m.idx[CH_POS][3*f+k] = (uint16_t)grid_corner(f, k, w);

    free(full_p); free(full_n); free(full_uv); free(full_c);
    return m;
}

void mesh_free(Mesh* m) {
    free(m->pos); free(m->nrm); free(m->uv); free(m->col);
    for (int ch = 0; ch < CH_N; ch++) free(m->idx[ch]);
    memset(m, 0, sizeof(*m));
}

void mesh_counts(const Mesh* m, int ch, int* counts, int via_vertex) {
    memset(counts, 0, (size_t)m->bank[ch] * sizeof(int));
    if (!via_vertex) {
        for (int k = 0; k < m->nf * 3; k++) {
            int s = m->idx[ch][k];
            if (s >= 0 && s < m->bank[ch]) counts[s]++;
        }
        return;
    }
    /* category error: count VERTEX ids from the position stream, report them
     * on rows 0..nv-1 of this channel's bank.  Rows >= nv stay zero. */
    for (int k = 0; k < m->nf * 3; k++) {
        int v = m->idx[CH_POS][k];
        if (v >= 0 && v < m->nv && v < m->bank[ch]) counts[v]++;
    }
}

Mesh mesh_prune(const Mesh* src, int* const live[CH_N]) {
    Mesh d; memset(&d, 0, sizeof(d));
    d.name = src->name; d.nv = src->nv; d.nf = src->nf;
    int* remap[CH_N];
    for (int ch = 0; ch < CH_N; ch++) {
        remap[ch] = malloc((size_t)src->bank[ch] * sizeof(int));
        int k = 0;
        for (int s = 0; s < src->bank[ch]; s++) remap[ch][s] = live[ch][s] ? k++ : 0;
        d.bank[ch] = k;
        d.live_n[ch] = k;
    }
    d.pos = calloc((size_t)(d.bank[CH_POS] + 1), sizeof(V3));
    d.nrm = calloc((size_t)(d.bank[CH_NRM] + 1), sizeof(V3));
    d.uv  = calloc((size_t)(d.bank[CH_UV] + 1) * 2, sizeof(float));
    d.col = calloc((size_t)(d.bank[CH_COL] + 1) * 3, sizeof(float));
    for (int s = 0; s < src->bank[CH_POS]; s++)
        if (live[CH_POS][s]) d.pos[remap[CH_POS][s]] = src->pos[s];
    for (int s = 0; s < src->bank[CH_NRM]; s++)
        if (live[CH_NRM][s]) d.nrm[remap[CH_NRM][s]] = src->nrm[s];
    for (int s = 0; s < src->bank[CH_UV]; s++)
        if (live[CH_UV][s]) {
            d.uv[2*remap[CH_UV][s]] = src->uv[2*s];
            d.uv[2*remap[CH_UV][s]+1] = src->uv[2*s+1];
        }
    for (int s = 0; s < src->bank[CH_COL]; s++)
        if (live[CH_COL][s]) {
            d.col[3*remap[CH_COL][s]+0] = src->col[3*s+0];
            d.col[3*remap[CH_COL][s]+1] = src->col[3*s+1];
            d.col[3*remap[CH_COL][s]+2] = src->col[3*s+2];
        }
    for (int ch = 0; ch < CH_N; ch++) {
        d.idx[ch] = malloc((size_t)src->nf * 3 * sizeof(uint16_t));
        for (int k = 0; k < src->nf * 3; k++) {
            int s = src->idx[ch][k];
            d.idx[ch][k] = (uint16_t)((s >= 0 && s < src->bank[ch]) ? remap[ch][s] : 0);
        }
        free(remap[ch]);
    }
    return d;
}

size_t mesh_bank_bytes(const Mesh* m, int ch) {
    return (size_t)m->bank[ch] * (size_t)CH_COMP[ch] * 4;
}
size_t mesh_attr_bytes(const Mesh* m) {
    size_t t = 0;
    for (int ch = 0; ch < CH_N; ch++) t += mesh_bank_bytes(m, ch);
    return t;
}
size_t mesh_index_bytes(const Mesh* m) { return (size_t)m->nf * 3 * 2 * CH_N; }

/* file layout: sample-style header, then the four index streams
 * (uint16, same width before and after prune), then the four banks.
 * The bank slot counts live in the header, so headers are 40 bytes both
 * before and after -- the index section width never changes. */
#define HDR 40

static void put_i32(unsigned char* p, int v) { memcpy(p, &v, 4); }
static int  get_i32(const unsigned char* p) { int v; memcpy(&v, p, 4); return v; }

int mesh_write(const Mesh* m, const char* path) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return -1;
    unsigned char hdr[HDR];
    memset(hdr, 0, sizeof(hdr));
    memcpy(hdr, "MMS2", 4);
    put_i32(hdr + 4, m->nv);
    put_i32(hdr + 8, m->nf);
    put_i32(hdr + 12, m->nf * 3);            /* stream length per channel */
    for (int ch = 0; ch < CH_N; ch++)
        put_i32(hdr + 16 + 4 * ch, m->bank[ch]);
    put_i32(hdr + 32, (int)mesh_attr_bytes(m));
    put_i32(hdr + 36, (int)mesh_index_bytes(m));
    fwrite(hdr, 1, sizeof(hdr), fp);
    for (int ch = 0; ch < CH_N; ch++)
        fwrite(m->idx[ch], sizeof(uint16_t), (size_t)m->nf * 3, fp);
    fwrite(m->pos, sizeof(V3), (size_t)m->bank[CH_POS], fp);
    fwrite(m->nrm, sizeof(V3), (size_t)m->bank[CH_NRM], fp);
    fwrite(m->uv, sizeof(float) * 2, (size_t)m->bank[CH_UV], fp);
    fwrite(m->col, sizeof(float) * 3, (size_t)m->bank[CH_COL], fp);
    fclose(fp);
    return 0;
}

Mesh mesh_read(const char* path) {
    Mesh m; memset(&m, 0, sizeof(m));
    FILE* fp = fopen(path, "rb");
    if (!fp) { m.nv = -1; return m; }
    unsigned char hdr[HDR];
    if (fread(hdr, 1, sizeof(hdr), fp) != sizeof(hdr)) { fclose(fp); m.nv = -1; return m; }
    if (memcmp(hdr, "MMS2", 4) != 0) { fclose(fp); m.nv = -1; return m; }
    m.nv = get_i32(hdr + 4);
    m.nf = get_i32(hdr + 8);
    int len = get_i32(hdr + 12);
    for (int ch = 0; ch < CH_N; ch++) { m.bank[ch] = get_i32(hdr + 16 + 4 * ch); m.live_n[ch] = m.bank[ch]; }
    m.pos = malloc((size_t)m.bank[CH_POS] * sizeof(V3));
    m.nrm = malloc((size_t)m.bank[CH_NRM] * sizeof(V3));
    m.uv  = malloc((size_t)m.bank[CH_UV] * 2 * sizeof(float));
    m.col = malloc((size_t)m.bank[CH_COL] * 3 * sizeof(float));
    for (int ch = 0; ch < CH_N; ch++) {
        m.idx[ch] = malloc((size_t)len * sizeof(uint16_t));
        if ((int)fread(m.idx[ch], sizeof(uint16_t), (size_t)len, fp) != len) {
            fclose(fp); m.nv = -1; return m;
        }
    }
    size_t got = 0;
    got += fread(m.pos, sizeof(V3), (size_t)m.bank[CH_POS], fp);
    got += fread(m.nrm, sizeof(V3), (size_t)m.bank[CH_NRM], fp);
    got += fread(m.uv, sizeof(float) * 2, (size_t)m.bank[CH_UV], fp);
    got += fread(m.col, sizeof(float) * 3, (size_t)m.bank[CH_COL], fp);
    if (got == 0) { fclose(fp); m.nv = -1; return m; }
    fclose(fp);
    return m;
}

/* ---------- exact check ---------- */

int attrs_equal(const Mesh* a, const Mesh* b) {
    if (a->nf != b->nf) return 0;
    for (int f = 0; f < a->nf; f++) {
        for (int k = 0; k < 3; k++) {
            int pa = a->idx[CH_POS][3*f+k], pb = b->idx[CH_POS][3*f+k];
            int na = a->idx[CH_NRM][3*f+k], nb = b->idx[CH_NRM][3*f+k];
            int ua = a->idx[CH_UV][3*f+k],  ub = b->idx[CH_UV][3*f+k];
            int ca = a->idx[CH_COL][3*f+k], cb = b->idx[CH_COL][3*f+k];
            if (memcmp(&a->pos[pa], &b->pos[pb], sizeof(V3)) != 0) return 0;
            if (memcmp(&a->nrm[na], &b->nrm[nb], sizeof(V3)) != 0) return 0;
            if (memcmp(&a->uv[2*ua], &b->uv[2*ub], 2 * sizeof(float)) != 0) return 0;
            if (memcmp(&a->col[3*ca], &b->col[3*cb], 3 * sizeof(float)) != 0) return 0;
        }
    }
    return 1;
}

/* ---------- render ---------- */

static V3 bary(float px, float py, const float* sx, const float* sy, int* ok) {
    V3 r = { 0, 0, 0 };
    float d = (sy[1]-sy[2])*(sx[0]-sx[2]) + (sx[2]-sx[1])*(sy[0]-sy[2]);
    if (fabsf(d) < 1e-12f) { *ok = 0; return r; }
    r.x = ((sy[1]-sy[2])*(px-sx[2]) + (sx[2]-sx[1])*(py-sy[2])) / d;
    r.y = ((sy[2]-sy[0])*(px-sx[2]) + (sx[0]-sx[2])*(py-sy[2])) / d;
    r.z = 1.0f - r.x - r.y;
    *ok = (r.x >= -1e-5f && r.y >= -1e-5f && r.z >= -1e-5f);
    return r;
}

void mesh_render(const Mesh* m, int W, int H, float* rgb) {
    size_t N = (size_t)W * H;
    float* depth = malloc(N * sizeof(float));
    for (size_t i = 0; i < N * 3; i++) rgb[i] = 0.06f;
    for (size_t i = 0; i < N; i++) depth[i] = -1e30f;
    V3 L = { 0.4f, 0.5f, 0.75f };
    float ll = sqrtf(L.x*L.x + L.y*L.y + L.z*L.z);
    L.x /= ll; L.y /= ll; L.z /= ll;
    for (int f = 0; f < m->nf; f++) {
        float sx[3], sy[3], sz[3];
        for (int k = 0; k < 3; k++) {
            V3 p = m->pos[m->idx[CH_POS][3*f+k]];
            float ca = 0.55f, sa = 0.55f;
            float X = p.x * ca + p.z * sa;
            float Z = -p.x * sa + p.z * ca;
            float Y = p.y;
            sx[k] = (X + 1.25f) / 2.5f * (float)(W - 1);
            sy[k] = (1.0f - (Y + 1.25f) / 2.5f) * (float)(H - 1);
            sz[k] = Z;
        }
        int minx = (int)floorf(fminf(sx[0], fminf(sx[1], sx[2])));
        int maxx = (int)ceilf (fmaxf(sx[0], fmaxf(sx[1], sx[2])));
        int miny = (int)floorf(fminf(sy[0], fminf(sy[1], sy[2])));
        int maxy = (int)ceilf (fmaxf(sy[0], fmaxf(sy[1], sy[2])));
        if (minx < 0) minx = 0;
        if (miny < 0) miny = 0;
        if (maxx >= W) maxx = W - 1;
        if (maxy >= H) maxy = H - 1;
        for (int y = miny; y <= maxy; y++) {
            for (int x = minx; x <= maxx; x++) {
                int ok = 0;
                V3 b = bary((float)x + 0.5f, (float)y + 0.5f, sx, sy, &ok);
                if (!ok) continue;
                float z = b.x*sz[0] + b.y*sz[1] + b.z*sz[2];
                if (z <= depth[y*W+x]) continue;
                depth[y*W+x] = z;
                float w[3] = { b.x, b.y, b.z };
                float cr = 0, cg = 0, cb = 0;
                for (int k = 0; k < 3; k++) {
                    int sn = m->idx[CH_NRM][3*f+k];
                    int sc = m->idx[CH_COL][3*f+k];
                    V3 n = m->nrm[sn];
                    float shade = 0.25f + 0.75f * fabsf(n.x*L.x + n.y*L.y + n.z*L.z);
                    cr += w[k] * m->col[3*sc+0] * shade;
                    cg += w[k] * m->col[3*sc+1] * shade;
                    cb += w[k] * m->col[3*sc+2] * shade;
                }
                size_t o = (size_t)(y*W+x) * 3;
                rgb[o+0] = cr; rgb[o+1] = cg; rgb[o+2] = cb;
            }
        }
    }
    free(depth);
}

int pixels_differ(const float* a, const float* b, int W, int H,
                  long* ndiff, float* max_delta) {
    long n = 0; float mx = 0;
    for (int i = 0; i < W * H; i++) {
        int d = 0;
        for (int c = 0; c < 3; c++) {
            float x = a[3*i+c], y = b[3*i+c];
            if (memcmp(&x, &y, sizeof(float)) != 0) d = 1;
            float dd = fabsf(x - y);
            if (dd > mx) mx = dd;
        }
        if (d) n++;
    }
    *ndiff = n; *max_delta = mx;
    return n > 0;
}

int write_ppm(const char* path, const float* rgb, int W, int H) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return -1;
    fprintf(fp, "P6\n%d %d\n255\n", W, H);
    for (int i = 0; i < W * H * 3; i++) {
        float v = rgb[i];
        if (v < 0) v = 0;
        if (v > 1) v = 1;
        unsigned char c = (unsigned char)(v * 255.0f + 0.5f);
        fwrite(&c, 1, 1, fp);
    }
    fclose(fp);
    return 0;
}

int write_ppm_diff(const char* path, const float* a, const float* b, int W, int H) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return -1;
    fprintf(fp, "P6\n%d %d\n255\n", W, H);
    for (int i = 0; i < W * H; i++) {
        float d = 0;
        for (int c = 0; c < 3; c++) {
            float dd = fabsf(a[3*i+c] - b[3*i+c]);
            if (dd > d) d = dd;
        }
        float v = d * 8.0f;
        if (v > 1) v = 1;
        unsigned char c = (unsigned char)(v * 255.0f + 0.5f);
        unsigned char z = 0;
        fwrite(&c, 1, 1, fp); fwrite(&z, 1, 1, fp); fwrite(&z, 1, 1, fp);
    }
    fclose(fp);
    return 0;
}
