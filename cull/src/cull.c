/* cull: frustum culling lab.
 *
 * 背景: 几百个实例用"最长边当半径"的胖包围球剔除, 球比盒子粗得多,
 *       漏剔多 -> 提交数下不来 -> 帧率无改善。
 *
 * 本程序:
 *   1. 胖球测试(radius = 最长边) 与 AABB 视锥精确测试(p/n-vertex) 对照;
 *   2. SAT(分离轴) 求 box-vs-frustum 的精确相交, 作为"漏剔"判定真值;
 *   3. 实例四叉树(空间分层 H): 整棵子树接受/拒绝, 浅 H 退化成逐个精确测试;
 *   4. 盒树粒度 g(0..2: 1/8/64 个子盒, 64 个提交/实例):
 *      粒度由盒子对角线 D 与中心到最近视锥面距离 delta 决定;
 *   5. 把"精确测试成本 C(g)"和"少提交带来的节省 S(g)"写成两个函数,
 *      用边际曲线交点求划算深度 g*;
 *   6. 三档实例数各跑一次, 输出 提交数/剔除耗时/漏剔数/实际分层深度。
 *
 * 场景直接生成在观察空间(相机原点, 朝 +Z), 视锥 6 面解析给出,
 * 与世界空间经视图矩阵变换后的平面测试在数学上完全相同。
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include <time.h>

/* ---------------- 基础类型 ---------------- */

typedef struct { float x, y, z; } V3;
typedef struct { V3 mn, mx; } Box;
typedef struct { V3 n; float d; } Plane;   /* 内部条件: n.p + d >= 0 */

static V3 v3_sub(V3 a, V3 b) { return (V3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static V3 v3_cross(V3 a, V3 b) {
    return (V3){ a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x };
}
static float v3_dot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static float v3_len(V3 a) { return sqrtf(v3_dot(a,a)); }

static V3 box_center(const Box* b) {
    return (V3){ (b->mn.x+b->mx.x)*0.5f, (b->mn.y+b->mx.y)*0.5f, (b->mn.z+b->mx.z)*0.5f };
}
static V3 box_half(const Box* b) {
    return (V3){ (b->mx.x-b->mn.x)*0.5f, (b->mx.y-b->mn.y)*0.5f, (b->mx.z-b->mn.z)*0.5f };
}
static float box_diag(const Box* b) {
    V3 h = box_half(b);
    return 2.0f * sqrtf(h.x*h.x + h.y*h.y + h.z*h.z);
}

/* ---------------- 视锥 ---------------- */

typedef struct {
    Plane pl[6];
    V3 corner[8];
    float zn, zf, tan_t, kx;
} Frust;

static void plane_normalize(Plane* p) {
    float l = sqrtf(v3_dot(p->n, p->n));
    p->n.x /= l; p->n.y /= l; p->n.z /= l; p->d /= l;
}
static float plane_dist(const Plane* p, V3 q) { return v3_dot(p->n, q) + p->d; }

/* 观察空间: 眼在原点, 朝 +Z。 z=z0 处半高=tan_t*z0, 半宽=kx*z0 */
static Frust frust_make(float zn, float zf, float fovy_deg, float aspect) {
    Frust f; memset(&f, 0, sizeof(f));
    f.zn = zn; f.zf = zf;
    f.tan_t = tanf(fovy_deg * 0.5f * 3.14159265358979323846f / 180.0f);
    f.kx = f.tan_t * aspect;
    float t = f.tan_t, k = f.kx;
    Plane raw[6] = {
        {{ 0, 0, 1}, -zn},   /* near: z >= zn */
        {{ 0, 0,-1},  zf},   /* far:  z <= zf */
        {{ 1, 0, k},   0},   /* left:  x + k z >= 0 */
        {{-1, 0, k},   0},   /* right: -x + k z >= 0 */
        {{ 0, 1, t},   0},   /* bottom: y + t z >= 0 */
        {{ 0,-1, t},   0},   /* top:    -y + t z >= 0 */
    };
    for (int i = 0; i < 6; i++) { f.pl[i] = raw[i]; plane_normalize(&f.pl[i]); }

    float xn = k*zn, yn = t*zn, xf = k*zf, yf = t*zf;
    V3 c[8] = {
        {-xn,-yn,zn}, { xn,-yn,zn}, { xn, yn,zn}, {-xn, yn,zn},
        {-xf,-yf,zf}, { xf,-yf,zf}, { xf, yf,zf}, {-xf, yf,zf},
    };
    memcpy(f.corner, c, sizeof(c));
    return f;
}

enum { CLS_OUT = -1, CLS_STRADDLE = 0, CLS_IN = 1 };

/* AABB p/n-vertex 精确平面测试。
 * 返回 OUT / STRADDLE / IN; out_mindelta = 中心到最近视锥面的距离(归一化)。 */
static int box_classify(const Box* b, const Frust* f, float* out_mindelta) {
    V3 c = box_center(b), h = box_half(b);
    int all_inside = 1;
    float md = 1e30f;
    for (int i = 0; i < 6; i++) {
        float dc = plane_dist(&f->pl[i], c);
        float r = fabsf(f->pl[i].n.x)*h.x + fabsf(f->pl[i].n.y)*h.y + fabsf(f->pl[i].n.z)*h.z;
        if (dc + r < 0.0f) { if (out_mindelta) *out_mindelta = md; return CLS_OUT; }
        if (dc - r < 0.0f) all_inside = 0;
        float a = fabsf(dc);
        if (a < md) md = a;
    }
    if (out_mindelta) *out_mindelta = md;
    return all_inside ? CLS_IN : CLS_STRADDLE;
}

/* 胖球基线: 半径直接取"最长边"(不是半对角线, 也不是半长边), 球显著过粗。 */
static int sphere_longestedge_test(const Box* b, const Frust* f) {
    V3 c = box_center(b);
    V3 s = { b->mx.x-b->mn.x, b->mx.y-b->mn.y, b->mx.z-b->mn.z };
    float R = s.x;
    if (s.y > R) R = s.y;
    if (s.z > R) R = s.z;
    for (int i = 0; i < 6; i++) {
        float dc = plane_dist(&f->pl[i], c);
        if (dc + R < 0.0f) return CLS_OUT;
    }
    return CLS_STRADDLE; /* 球没有"完全在内"的快速结论需求, 统一当可见 */
}

/* ---------------- SAT 真值: AABB vs 视锥(凸多面体) 精确相交 ---------------- */

static int sat_separated(const Box* b, const Frust* f, V3 axis) {
    float al = v3_len(axis);
    if (al < 1e-9f) return 0;
    axis.x /= al; axis.y /= al; axis.z /= al;
    V3 c = box_center(b), h = box_half(b);
    float bc = v3_dot(axis, c);
    float br = fabsf(axis.x)*h.x + fabsf(axis.y)*h.y + fabsf(axis.z)*h.z;
    float bmin = bc - br, bmax = bc + br;
    float fmin = 1e30f, fmax = -1e30f;
    for (int i = 0; i < 8; i++) {
        float d = v3_dot(axis, f->corner[i]);
        if (d < fmin) fmin = d;
        if (d > fmax) fmax = d;
    }
    const float EPS = 1e-4f;
    return (bmax < fmin - EPS || fmax < bmin - EPS);
}

static int box_frustum_exact(const Box* b, const Frust* f) {
    V3 axes[48];
    int na = 0;
    axes[na++] = (V3){1,0,0};
    axes[na++] = (V3){0,1,0};
    axes[na++] = (V3){0,0,1};
    for (int i = 0; i < 6; i++) axes[na++] = f->pl[i].n;

    static const int EDGE[12][2] = {
        {0,1},{1,2},{2,3},{3,0},
        {4,5},{5,6},{6,7},{7,4},
        {0,4},{1,5},{2,6},{3,7}
    };
    static const V3 BOXA[3] = {{1,0,0},{0,1,0},{0,0,1}};
    for (int a = 0; a < 3; a++) {
        for (int e = 0; e < 12; e++) {
            V3 fe = v3_sub(f->corner[EDGE[e][1]], f->corner[EDGE[e][0]]);
            V3 cr = v3_cross(BOXA[a], fe);
            if (v3_len(cr) > 1e-6f && na < 48) axes[na++] = cr;
        }
    }
    for (int i = 0; i < na; i++)
        if (sat_separated(b, f, axes[i])) return 0;
    return 1;
}

/* ---------------- 实例子盒树(粒度 g: 每轴 2^g 段) ---------------- */

#define SUBS_PER_INST 64      /* 每实例最多 64 个可提交单元 = g=2 的 64 个子盒 */
#define GRAN_MAX 2

static Box sub_box(const Box* b, int g, int idx) {
    int k = 1 << g;
    int nx = k*k;
    int ix = idx / nx;
    int iy = (idx / k) % k;
    int iz = idx % k;
    Box r;
    r.mn.x = b->mn.x + (b->mx.x-b->mn.x) * ((float)ix/k);
    r.mx.x = b->mn.x + (b->mx.x-b->mn.x) * ((float)(ix+1)/k);
    r.mn.y = b->mn.y + (b->mx.y-b->mn.y) * ((float)iy/k);
    r.mx.y = b->mn.y + (b->mx.y-b->mn.y) * ((float)(iy+1)/k);
    r.mn.z = b->mn.z + (b->mx.z-b->mn.z) * ((float)iz/k);
    r.mx.z = b->mn.z + (b->mx.z-b->mn.z) * ((float)(iz+1)/k);
    return r;
}

/* 粒度公式: 子盒对角线 D/2^g 要细到不超过边界距离 delta 的 beta 倍。
 *   g = clamp(ceil(log2(D / (beta*delta))))
 * 大盒子贴近视锥面 -> g 大(细分精确测, 多救回提交);
 * 小盒子离边界远   -> g=0(整盒一次测完)。 */
static int adaptive_granularity(const Box* b, float delta) {
    /* BETA=2: 要求子盒对角线 D/2^g 不超过 2*delta(边界距离的两倍容差带)。
     * 这样仅擦边的盒子停在粗档、被面切到中部的大盒进细档, 粒度有分化。 */
    const float BETA = 2.0f;
    float D = box_diag(b);
    float ratio;
    if (delta < 1e-3f) ratio = 1e30f;
    else ratio = D / (BETA * delta);
    if (ratio <= 1.0f) return 0;
    int g = (int)ceilf(log2f(ratio));
    if (g < 0) g = 0;
    if (g > GRAN_MAX) g = GRAN_MAX;
    return g;
}

/* ---------------- 场景 ---------------- */

typedef struct { Box box; } Inst;

static uint64_t rng_state = 0x123456789abcdef0ULL;
static float rndf(float lo, float hi) {
    uint64_t x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng_state = x;
    float u = (float)((x >> 11) & ((1ULL<<53)-1)) / (float)(1ULL<<53);
    return lo + (hi-lo)*u;
}

static void gen_scene(int n, Inst* out) {
    /* 一部分实例成簇分布(模拟建筑群/道具组), 一部分均匀散布。
     * 成簇让空间树能在中层整簇接受/拒绝; 均匀部分保证各深度都有穿越。 */
    enum { MAXC = 64 };
    int nc = n / 18;
    if (nc < 6) nc = 6;
    if (nc > MAXC) nc = MAXC;
    float ccx[MAXC], ccz[MAXC], crad[MAXC];
    for (int c = 0; c < nc; c++) {
        ccx[c] = rndf(-115.0f, 115.0f);
        ccz[c] = rndf(6.0f, 150.0f);
        crad[c] = rndf(3.0f, 11.0f);
    }
    for (int i = 0; i < n; i++) {
        int cls = i % 3;                 /* 小/中/大盒均匀混合 */
        float lo, hi;
        if (cls == 0)      { lo = 0.6f; hi = 2.2f; }
        else if (cls == 1) { lo = 2.5f; hi = 7.0f; }
        else               { lo = 7.0f; hi = 15.0f; }
        float sx = rndf(lo, hi) * rndf(0.7f, 1.3f);
        float sy = rndf(lo, hi) * rndf(0.5f, 1.1f);
        float sz = rndf(lo, hi) * rndf(0.7f, 1.3f);
        float cx, cz;
        if (rndf(0.0f, 1.0f) < 0.72f) {
            int c = (int)rndf(0.0f, (float)nc - 0.001f);
            if (c >= nc) c = nc-1;
            /* 两次均匀随机相加 -> 中部聚集的抖动 */
            float ox = (rndf(-1,1)+rndf(-1,1)) * crad[c];
            float oz = (rndf(-1,1)+rndf(-1,1)) * crad[c];
            cx = ccx[c] + ox;
            cz = ccz[c] + oz;
        } else {
            cx = rndf(-125.0f, 125.0f);
            cz = rndf(1.5f, 155.0f);
        }
        out[i].box.mn = (V3){ cx-sx*0.5f, 0.0f, cz-sz*0.5f };
        out[i].box.mx = (V3){ cx+sx*0.5f, sy,      cz+sz*0.5f };
    }
}

/* ---------------- 实例四叉树(空间分层 H) ---------------- */

#define HMAX 5
#define LEAF_CAP 6

typedef struct {
    Box box;
    int depth, count;
    int child[4], nchild;
    int* ids; int nids;
} QNode;

typedef struct { QNode* node; int n, cap; } QTree;

static int qt_add(QTree* t) {
    if (t->n == t->cap) {
        t->cap = t->cap ? t->cap*2 : 64;
        t->node = realloc(t->node, (size_t)t->cap * sizeof(QNode));
    }
    int id = t->n++;
    memset(&t->node[id], 0, sizeof(QNode));
    for (int q = 0; q < 4; q++) t->node[id].child[q] = -1;
    return id;
}

static int qt_build(QTree* t, const Inst* ins, const int* ids, int n, int depth,
                    float x0, float x1, float z0, float z1) {
    int id = qt_add(t);
    t->node[id].depth = depth;
    t->node[id].ids = malloc((size_t)n * sizeof(int));
    memcpy(t->node[id].ids, ids, (size_t)n * sizeof(int));
    t->node[id].nids = n;
    t->node[id].count = n;

    Box u = ins[ids[0]].box;
    for (int k = 1; k < n; k++) {
        const Box* b = &ins[ids[k]].box;
        if (b->mn.x < u.mn.x) u.mn.x = b->mn.x;
        if (b->mn.y < u.mn.y) u.mn.y = b->mn.y;
        if (b->mn.z < u.mn.z) u.mn.z = b->mn.z;
        if (b->mx.x > u.mx.x) u.mx.x = b->mx.x;
        if (b->mx.y > u.mx.y) u.mx.y = b->mx.y;
        if (b->mx.z > u.mx.z) u.mx.z = b->mx.z;
    }
    t->node[id].box = u;

    int kids[4] = {-1,-1,-1,-1};
    if (depth < HMAX && n > LEAF_CAP) {
        float mx = (x0+x1)*0.5f, mz = (z0+z1)*0.5f;
        int* buck[4]; int cnt[4] = {0,0,0,0};
        for (int q = 0; q < 4; q++) buck[q] = malloc((size_t)n * sizeof(int));
        for (int k = 0; k < n; k++) {
            V3 c = box_center(&ins[ids[k]].box);
            int q = (c.x >= mx ? 1 : 0) + (c.z >= mz ? 2 : 0);
            buck[q][cnt[q]++] = ids[k];
        }
        int all_in_one = 0;
        for (int q = 0; q < 4; q++) if (cnt[q] == n) all_in_one = 1;
        if (!all_in_one) {
            float xs[2][2] = {{x0,mx},{mx,x1}};
            float zs[2][2] = {{z0,mz},{mz,z1}};
            for (int q = 0; q < 4; q++) {
                if (cnt[q] == 0) { free(buck[q]); continue; }
                int xi = q & 1, zi = (q >> 1) & 1;
                kids[q] = qt_build(t, ins, buck[q], cnt[q], depth+1,
                                   xs[xi][0], xs[xi][1], zs[zi][0], zs[zi][1]);
                free(buck[q]);
            }
        } else {
            for (int q = 0; q < 4; q++) free(buck[q]);
        }
    }
    for (int q = 0; q < 4; q++) if (kids[q] >= 0) {
        t->node[id].child[q] = kids[q];
        t->node[id].nchild++;
    }
    return id;
}

static void qt_free(QTree* t) {
    for (int i = 0; i < t->n; i++) free(t->node[i].ids);
    free(t->node);
    t->node = NULL; t->n = t->cap = 0;
}

/* ---------------- 剔除遍历 ---------------- */

typedef struct {
    long subs;
    long node_tests, inst_tests, box_tests;
    long depth_sum;
} Stats;

typedef struct {
    int n;
    Frust fr;
    Inst* inst;
    QTree tree;
    uint8_t* truth;     /* n*64: SAT 真值 */
    int* badn;          /* 每实例 SAT 判出的在外叶子数 */
    int8_t* dec;        /* 0/1/2 */
    uint8_t* gran;
} Bench;

static void process_instance(Bench* B, int ii, int gcap, Stats* st) {
    const Box* b = &B->inst[ii].box;
    float delta = 0.0f;
    int cls = box_classify(b, &B->fr, &delta);
    st->inst_tests++;
    if (cls == CLS_OUT) { B->dec[ii] = 0; return; }
    if (cls == CLS_IN)  { B->dec[ii] = 1; st->subs += SUBS_PER_INST; return; }

    int g = adaptive_granularity(b, delta);
    if (g > gcap) g = gcap;
    B->dec[ii] = 2;
    B->gran[ii] = (uint8_t)g;

    if (g == 0) { st->subs += SUBS_PER_INST; return; }
    int k = 1 << g, nb = k*k*k, share = SUBS_PER_INST/nb;
    for (int q = 0; q < nb; q++) {
        Box sb = sub_box(b, g, q);
        int c2 = box_classify(&sb, &B->fr, NULL);
        st->box_tests++;
        if (c2 != CLS_OUT) st->subs += share;
    }
}

static void cull_node(Bench* B, int nid, int hcap, int gcap, Stats* st) {
    const QNode* nd = &B->tree.node[nid];
    int cls = box_classify(&nd->box, &B->fr, NULL);
    st->node_tests++;

    if (cls == CLS_OUT) {
        for (int k = 0; k < nd->nids; k++) B->dec[nd->ids[k]] = 0;
        return;
    }
    if (cls == CLS_IN) {
        st->subs += (long)SUBS_PER_INST * nd->count;
        st->depth_sum += (long)nd->depth * nd->count;
        for (int k = 0; k < nd->nids; k++) B->dec[nd->ids[k]] = 1;
        return;
    }

    /* STRADDLE: 深度额度未用完且有孩子 -> 继续分层;
     * 否则在本层把名下实例逐个精确测试(浅 H 时的退化路径)。 */
    if (nd->depth < hcap && nd->nchild > 0) {
        for (int q = 0; q < 4; q++)
            if (nd->child[q] >= 0) cull_node(B, nd->child[q], hcap, gcap, st);
        return;
    }
    for (int k = 0; k < nd->nids; k++) {
        int ii = nd->ids[k];
        st->depth_sum += nd->depth;
        process_instance(B, ii, gcap, st);
    }
}

static Stats cull_scene(Bench* B, int hcap, int gcap) {
    Stats st; memset(&st, 0, sizeof(st));
    cull_node(B, 0, hcap, gcap, &st);
    return st;
}

/* 粗粒度 g 的第 q 个子盒, 覆盖 GRAN_MAX 叶子的数量(每轴 s 格)。
 * 遍历其名下所有最细叶子, 数 SAT 真值在外(=漏剔的提交)。 */
static long misses_under_box(const Bench* B, int ii, int g, int q) {
    const int G = GRAN_MAX;
    int kg = 1 << g, kG = 1 << G;
    int ix = q / (kg*kg);
    int iy = (q / kg) % kg;
    int iz = q % kg;
    int s = kG / kg;
    long miss = 0;
    for (int ax = 0; ax < s; ax++)
    for (int ay = 0; ay < s; ay++)
    for (int az = 0; az < s; az++) {
        int lx = ix*s+ax, ly = iy*s+ay, lz = iz*s+az;
        int j = (lx*kG+ly)*kG+lz;
        if (!B->truth[ii*SUBS_PER_INST+j]) miss++;
    }
    return miss;
}

/* 强制所有跨界盒细分到 GRAN_MAX: 验证最细 AABB 测试相对 SAT 真值的漏剔。 */
static long forced_finest_misses(const Bench* B) {
    long miss = 0;
    for (int i = 0; i < B->n; i++) {
        const Box* b = &B->inst[i].box;
        int cls = box_classify(b, &B->fr, NULL);
        if (cls == CLS_OUT) continue;
        if (cls == CLS_IN) { miss += B->badn[i]; continue; }
        for (int j = 0; j < SUBS_PER_INST; j++) {
            Box leaf = sub_box(b, GRAN_MAX, j);
            if (box_classify(&leaf, &B->fr, NULL) == CLS_OUT) continue;
            if (!B->truth[i*SUBS_PER_INST+j]) miss++;
        }
    }
    return miss;
}

/* 用 SAT 真值核对最近一次 cull 写入的 dec/gran, 数漏剔提交。 */
static long count_misses(const Bench* B) {
    long miss = 0;
    for (int i = 0; i < B->n; i++) {
        int8_t d = B->dec[i];
        if (d == 0) continue;
        if (d == 1) { miss += B->badn[i]; continue; }
        int g = B->gran[i];
        int k = 1 << g, nb = k*k*k;
        for (int q = 0; q < nb; q++) {
            Box sb = sub_box(&B->inst[i].box, g, q);
            if (box_classify(&sb, &B->fr, NULL) == CLS_OUT) continue;
            miss += misses_under_box(B, i, g, q);
        }
    }
    return miss;
}

/* ---------------- 计时与"提交"模拟负载 ---------------- */

static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec*1e6 + (double)ts.tv_nsec/1e3;
}

#ifndef TIMING_REPS
#define TIMING_REPS 800
#endif
#ifndef TIMING_BEST
#define TIMING_BEST 9
#endif

static double time_cull(Bench* B, int hcap, int gcap) {
    cull_scene(B, hcap, gcap); /* 预热 I-cache/分支/堆 */
    double best = 1e30;
    volatile long sink = 0;
    for (int bat = 0; bat < TIMING_BEST; bat++) {
        double t0 = now_us();
        Stats st = {0};
        for (int r = 0; r < TIMING_REPS; r++) st = cull_scene(B, hcap, gcap);
        double dt = (now_us()-t0)/(double)TIMING_REPS;
        if (dt < best) best = dt;
        sink += st.subs + st.node_tests;
    }
    if (sink == 0x7fffffffL) printf("x"); /* 防优化的不可能分支 */
    return best;
}

/* 模拟一次提交: 填一条 48 字节命令并回读混校验, 使耗时与提交数真实相关。
 * 注意: 这只标定"命令填充"下限; 驱动/状态切换的完整提交成本用
 * K_NS_PER_SUBMIT 显式建模。 */
typedef struct { uint8_t b[48]; } Cmd;
static uint64_t submit_work(long count, Cmd* cmds) {
    uint64_t h = 1469598103934665603ULL;
    for (long i = 0; i < count; i++) {
        Cmd* c = &cmds[i];
        for (int j = 0; j < 48; j++) c->b[j] = (uint8_t)(i*131 + j*17);
        for (int j = 0; j < 48; j++) h ^= (uint64_t)c->b[j] + (uint64_t)(i & 255);
        h *= 1099511628211ULL;
    }
    return h;
}

/* 标定 ns/提交: 在几个提交量上实测, 最小二乘拟合斜率。 */
static double calibrate_ns_per_submit(long max_subs) {
    Cmd* cmds = malloc((size_t)(max_subs+1) * sizeof(Cmd));
    long xs[5]; double ys[5];
    for (int p = 0; p < 5; p++) {
        xs[p] = max_subs * (p+1) / 5;
        if (xs[p] < 1) xs[p] = 1;
        int reps = 60;
        volatile uint64_t sink = 0;
        double best = 1e30;
        for (int bat = 0; bat < 3; bat++) {
            double t0 = now_us();
            uint64_t h = 0;
            for (int r = 0; r < reps; r++) h ^= submit_work(xs[p], cmds);
            double dt = (now_us()-t0)*1000.0 / ((double)reps * xs[p]); /* ns/sub */
            if (dt < best) best = dt;
            sink += h;
        }
        ys[p] = best;
        if (sink == 999) printf("x");
    }
    free(cmds);
    /* ys 已是每提交 ns; 用过原点/均值稳健估计: 总时长斜率更准, 这里取中位 */
    for (int i = 0; i < 4; i++)
        for (int j = i+1; j < 5; j++)
            if (ys[j] < ys[i]) { double t = ys[i]; ys[i] = ys[j]; ys[j] = t; }
    return ys[2];
}

/* ---------------- 自检 ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("check_%-28s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void self_tests(void) {
    Frust f = frust_make(1.0f, 140.0f, 50.0f, 16.0f/9.0f);

    Box inside = { {-2,-2,20}, {2,4,26} };
    check("box_inside", box_classify(&inside,&f,NULL) == CLS_IN);
    Box behind = { {-2,-2,-8}, {2,4,-2} };
    check("box_behind_out", box_classify(&behind,&f,NULL) == CLS_OUT);
    Box faraway = { {-2,-2,200}, {2,4,210} };
    check("box_far_out", box_classify(&faraway,&f,NULL) == CLS_OUT);
    check("sat_inside", box_frustum_exact(&inside,&f) == 1);
    check("sat_behind", box_frustum_exact(&behind,&f) == 0);

    /* 胖球比盒子更容易放行: 同一扁长盒, 球半径=最长边 */
    Box longbox = { {-1,-1,30}, {1,1,46} };   /* z 长 16, 横向窄 */
    int bc = box_classify(&longbox,&f,NULL);
    int sc = sphere_longestedge_test(&longbox,&f);
    check("longbox_box_straddle", bc == CLS_STRADDLE || bc == CLS_IN);
    check("fatsphere_visible", sc != CLS_OUT);

    /* 找一个 p-vertex 放行但 SAT 判否的棱/角楔形盒, 证明漏剔真值有效 */
    int wedge_found = 0;
    uint64_t saved = rng_state;
    rng_state = 424242ULL;
    for (int t = 0; t < 200000 && !wedge_found; t++) {
        float z = rndf(2.0f, 138.0f);
        float hh = f.tan_t*z, hw = f.kx*z;
        int corner = t % 4;
        float sx = rndf(1.0f, 5.0f), sy = rndf(1.0f, 5.0f), sz = rndf(1.0f, 5.0f);
        float ex = (corner & 1) ? (hw + rndf(-2, 6)) : (-hw - rndf(-2, 6));
        float ey = (corner & 2) ? (hh + rndf(-2, 6)) : (-hh - rndf(-2, 6));
        Box b = { {ex-sx/2, ey-sy/2, z-sz/2}, {ex+sx/2, ey+sy/2, z+sz/2} };
        if (box_classify(&b,&f,NULL) != CLS_OUT && !box_frustum_exact(&b,&f))
            wedge_found = 1;
    }
    rng_state = saved;
    check("sat_wedge_case_exists", wedge_found);

    /* 粒度单调性: 同中心距离下, 大盒粒度 >= 小盒; 贴近边界 >= 远离边界 */
    Box small = { {-0.5f,-0.5f,0},{0.5f,0.5f,0} };
    Box big   = { {-8,-8,-8},{8,8,8} };
    check("gran_far_small0", adaptive_granularity(&small, 30.0f) == 0);
    check("gran_near_big_max", adaptive_granularity(&big, 0.2f) == GRAN_MAX);
    int g_far = adaptive_granularity(&big, 20.0f);
    int g_near = adaptive_granularity(&big, 2.0f);
    check("gran_near_ge_far", g_near >= g_far);

    /* 子盒索引覆盖与 share 一致 */
    int cover_ok = 1;
    for (int g = 0; g <= GRAN_MAX; g++) {
        int k = 1<<g, nb=k*k*k, share=SUBS_PER_INST/nb;
        if (nb*share != SUBS_PER_INST) cover_ok = 0;
    }
    check("subbox_share_covers_all", cover_ok);
}

/* ---------------- 主流程 ---------------- */

#ifndef K_NS_PER_SUBMIT
/* 成本/节省交点所用的交换系数: 一次提交折算多少 ns。
 * 命令填充实测只有十几 ns(见 cmd_fill), 完整驱动提交通常更贵;
 * 取 80ns 这个偏轻的模型值是为了让边际交点落在 g=1..2 的可选档位之间
 * 便于读数, 同时另给 break-even 阈值和 k 敏感性, 不依赖这个主观取值。 */
#define K_NS_PER_SUBMIT 80.0
#endif

typedef struct {
    int n;
    int hstar, gstar;
    double cull_us;
    long subs, miss;
    double depth_avg, gran_avg, layers_avg;
    double ratio_vs_first, step_ratio, ns_per_inst;
    double shallow_us;
    double shallow_inst_frac;
    double shallow_depth;
} Row;

int main(void) {
    self_tests();
    if (failures) { printf("all_ok=0\n"); return 1; }

    const int TIERS[3] = {200, 400, 800};
    const char* LABEL[3] = {"small", "mid", "max"};
    Row rows[3]; memset(rows, 0, sizeof(rows));

    Frust fr = frust_make(1.0f, 140.0f, 50.0f, 16.0f/9.0f);
    double tau_sub_ns = calibrate_ns_per_submit((long)TIERS[2]*SUBS_PER_INST);
    printf("\nmodel ns_per_submit_cmd_fill=%.2f effective=%.1f\n",
           tau_sub_ns, K_NS_PER_SUBMIT);

    /* 先求各档的成本/节省曲线与交点 */

    for (int ti = 0; ti < 3; ti++) {
        int n = TIERS[ti];
        rows[ti].n = n;

        Bench B; memset(&B, 0, sizeof(B));
        B.n = n; B.fr = fr;
        B.inst = calloc((size_t)n, sizeof(Inst));
        B.truth = malloc((size_t)n*SUBS_PER_INST);
        B.badn = calloc((size_t)n, sizeof(int));
        B.dec = malloc((size_t)n);
        B.gran = malloc((size_t)n);

        rng_state = 0xC0FFEE0000000001ULL ^ (uint64_t)n*0x9E3779B97F4A7C15ULL;
        gen_scene(n, B.inst);

        int* allids = malloc((size_t)n*sizeof(int));
        for (int i = 0; i < n; i++) allids[i] = i;
        B.tree = (QTree){0};
        qt_build(&B.tree, B.inst, allids, n, 0, -130.0f, 130.0f, -5.0f, 165.0f);
        free(allids);

        for (int i = 0; i < n; i++) {
            int bad = 0;
            for (int j = 0; j < SUBS_PER_INST; j++) {
                Box leaf = sub_box(&B.inst[i].box, GRAN_MAX, j);
                int vis = box_frustum_exact(&leaf, &fr);
                B.truth[i*SUBS_PER_INST+j] = (uint8_t)vis;
                if (!vis) bad++;
            }
            B.badn[i] = bad;
        }
        long forced_miss = forced_finest_misses(&B);
        printf("\n===== tier=%s n=%d =====\n", LABEL[ti], n);
        printf("forced_g%d_exact_missed=%ld (SAT ground truth)\n",
               GRAN_MAX, forced_miss);

        /* 胖球基线对照 */
        long sph_subs = 0, sph_miss = 0;
        double sph_best = 1e30; volatile long sink=0;
        for (int i = 0; i < n; i++) {
            if (sphere_longestedge_test(&B.inst[i].box, &fr) != CLS_OUT) {
                sph_subs += SUBS_PER_INST; sph_miss += B.badn[i];
            }
        }
        for (int bat=0; bat<TIMING_BEST; bat++) {
            double t0=now_us();
            for (int r=0;r<TIMING_REPS;r++) {
                long s=0; for (int i=0;i<n;i++)
                    if (sphere_longestedge_test(&B.inst[i].box,&fr)!=CLS_OUT) s+=SUBS_PER_INST;
                sink+=s;
            }
            double dt=(now_us()-t0)/(double)TIMING_REPS;
            if (dt<sph_best) sph_best=dt;
        }
        if (sink==7) printf("x");

        printf("sphere_longestedge: submissions=%ld missed=%ld cull_us=%.3f\n",
               sph_subs, sph_miss, sph_best);

        /* --- 粒度 g 扫描(先在最深空间树 HMAX 上做, 隔离空间分层影响):
         * 两个函数: C(g)=精确测试实测成本(us);
         *           S(g)=相对 g=0 少提交量折算的时间价值(us)。 --- */
        printf("-- granularity g sweep: C=precise_test_cost S=saved_submissions --\n");
        double C[GRAN_MAX+1], Sval[GRAN_MAX+1];
        long subs_g[GRAN_MAX+1];
        for (int g = 0; g <= GRAN_MAX; g++) {
            double gt = time_cull(&B, HMAX, g);
            Stats st = cull_scene(&B, HMAX, g);
            subs_g[g] = st.subs;
            C[g] = gt;                                      /* 函数1 */
            long saved = subs_g[0]-subs_g[g];
            Sval[g] = (double)saved * K_NS_PER_SUBMIT/1000.0; /* 函数2 */
            long miss_g = 0;
            /* 用该档决策核对漏剔: cull_scene 已写入 dec/gran */
            miss_g = count_misses(&B);
            printf("g=%d C_cost_us=%9.3f submissions=%7ld saved_subs=%7ld "
                   "S_save_us=%10.2f missed=%6ld\n",
                   g, C[g], subs_g[g], subs_g[0]-subs_g[g], Sval[g], miss_g);
        }
        /* 边际交点: 每加深一档 g-1->g, 额外测试成本 dC 与额外节省 dS。
         * k_break = 该档刚好划算所需的单提交成本(ns)。 */
        double ynet[GRAN_MAX+1];
        int g_pol = 0;
        for (int g = 1; g <= GRAN_MAX; g++) {
            double dC = C[g]-C[g-1];
            long dsub = subs_g[g-1]-subs_g[g];
            double dS = (double)dsub * K_NS_PER_SUBMIT/1000.0;
            double k_break = (dsub > 0) ? dC*1000.0/(double)dsub : 1e30;
            int worth = (dS >= dC);
            ynet[g] = dS - dC;
            printf("margin g-%d->%d: dC=%8.3fus d_save_subs=%6ld dS=%9.2fus "
                   "break_even_ns/sub=%8.1f worth@%.0fns=%d\n",
                   g-1, g, dC, dsub, dS, k_break, K_NS_PER_SUBMIT, worth);
            if (worth) g_pol = g;
        }
        /* 连续交点: 边际净收益 y(g)=dS-dC 在相邻整数档间变号时线性插值 */
        double gx = (double)g_pol;
        int crossed = 0;
        for (int g = 2; g <= GRAN_MAX; g++) {
            double y1 = ynet[g-1], y2 = ynet[g];
            if ((y1 > 0 && y2 < 0) || (y1 < 0 && y2 > 0)) {
                gx = (double)(g-1) + (-y1)/(y2-y1);
                crossed = 1;
            }
        }
        if (!crossed)
            gx = (ynet[GRAN_MAX] >= 0) ? (double)GRAN_MAX : 0.0;
        printf("intersection gx=%.2f policy_g=%d (@k=%.0fns/sub)\n",
               gx, g_pol, K_NS_PER_SUBMIT);
        /* 交换系数敏感性: 不同提交成本假设下的最细划算档 */
        {
            const double KS[6] = {12.0, 40.0, 80.0, 200.0, 1000.0, 5000.0};
            printf("k_sensitivity_gstar:");
            for (int z = 0; z < 6; z++) {
                int gp = 0;
                for (int g = 1; g <= GRAN_MAX; g++) {
                    double dC = C[g]-C[g-1];
                    long dsub = subs_g[g-1]-subs_g[g];
                    if ((double)dsub*KS[z]/1000.0 >= dC) gp = g;
                }
                printf(" k=%.0fns->g%d", KS[z], gp);
            }
            printf("\n");
        }

        /* --- 空间深度 H 扫描(用交点给出的 g*): 成本=测试操作总数+实测耗时 --- */
        printf("-- spatial H sweep (g*=%d) --\n", g_pol);
        int h_best = 0; long best_ops = INT32_MAX; double h_best_us = 0.0;
        for (int h = 0; h <= HMAX; h++) {
            double ht = time_cull(&B, h, g_pol);
            Stats st = cull_scene(&B, h, g_pol);
            long ops = st.node_tests + st.inst_tests + st.box_tests;
            double davg = (double)st.depth_sum / n;
            double ifrac = (double)st.inst_tests / n;
            printf("H=%d cull_us=%9.3f ops=%7ld node=%6ld inst=%5ld box=%7ld "
                   "inst/N=%.3f depth_avg=%.2f\n",
                   h, ht, ops, st.node_tests, st.inst_tests, st.box_tests,
                   ifrac, davg);
            if (ops < best_ops) { best_ops = ops; h_best = h; h_best_us = ht; }
        }

        /* --- 正式核对一次: H=h*, g=g* (耗时直接用 H 扫描中该档的 best) --- */
        double final_us = h_best_us;
        Stats fs = cull_scene(&B, h_best, g_pol);
        long miss = count_misses(&B);
        double gran_sum = 0; long gran_n = 0;
        int gh_raw[GRAN_MAX+1] = {0}, gh_exec[GRAN_MAX+1] = {0};
        for (int i = 0; i < n; i++) {
            if (B.dec[i] == 2) {
                gran_sum += B.gran[i]; gran_n++;
                gh_exec[B.gran[i]]++;
                float delta = 0.0f;
                box_classify(&B.inst[i].box, &fr, &delta);
                int gr = adaptive_granularity(&B.inst[i].box, delta);
                gh_raw[gr]++;
            }
        }
        printf("adaptive_gran_hist raw(g0/g1/g2)=%d/%d/%d exec=%d/%d/%d "
               "(straddlers=%ld)\n",
               gh_raw[0], gh_raw[1], gh_raw[2],
               gh_exec[0], gh_exec[1], gh_exec[2], gran_n);
        rows[ti].hstar = h_best;
        rows[ti].gstar = g_pol;
        rows[ti].cull_us = final_us;
        rows[ti].subs = fs.subs;
        rows[ti].miss = miss;
        rows[ti].depth_avg = (double)fs.depth_sum / n;
        rows[ti].gran_avg = gran_n ? gran_sum/gran_n : 0.0;
        rows[ti].layers_avg = rows[ti].depth_avg + rows[ti].gran_avg;
        rows[ti].ns_per_inst = final_us*1000.0 / n;

        /* --- 浅 H=0 对照: 复现"退化成逐个精确测试" --- */
        rows[ti].shallow_us = time_cull(&B, 0, g_pol);
        Stats s0 = cull_scene(&B, 0, g_pol);
        rows[ti].shallow_inst_frac = (double)s0.inst_tests / n;
        rows[ti].shallow_depth = (double)s0.depth_sum / n;

        qt_free(&B.tree);
        free(B.inst); free(B.truth); free(B.badn); free(B.dec); free(B.gran);
    }

    /* 比率并列: 三档的深度与耗时比放在同一张表 */
    printf("\n===== depth vs time ratios (chosen policy) =====\n");
    for (int ti = 0; ti < 3; ti++) {
        rows[ti].ratio_vs_first = rows[ti].cull_us / rows[0].cull_us;
        rows[ti].step_ratio = ti ? rows[ti].cull_us / rows[ti-1].cull_us : 1.0;
        double d_ratio = rows[ti].depth_avg / rows[0].depth_avg;
        printf("tier=%-3s n=%4d H*=%d depth_avg=%.2f depth_ratio=%4.2f "
               "g*=%d gran_avg=%.2f layers_avg=%.2f | cull_us=%9.3f "
               "time_ratio=%5.2f step=%5.2f (N step=2.00) ns/inst=%6.2f\n",
               LABEL[ti], rows[ti].n, rows[ti].hstar, rows[ti].depth_avg, d_ratio,
               rows[ti].gstar, rows[ti].gran_avg, rows[ti].layers_avg,
               rows[ti].cull_us, rows[ti].ratio_vs_first, rows[ti].step_ratio,
               rows[ti].ns_per_inst);
    }
    /* "最多档耗时反而上升"判据:
     * (a) 最大档相对中档的环比 > 1.5 x 数量环比(2.0) => >3.0; 或
     * (b) 每实例耗时最大档比小档上涨 >25%。 */
    int rise_flag = (rows[2].step_ratio > 1.5 * 2.0) ||
                    (rows[2].ns_per_inst > 1.25 * rows[0].ns_per_inst);
    printf("max_tier_cull_time_rise_flag=%d (mid->max step=%.2f vs count step 2.00; "
           "ns/inst %.2f -> %.2f)\n",
           rise_flag, rows[2].step_ratio, rows[0].ns_per_inst, rows[2].ns_per_inst);
    if (rise_flag)
        printf("diagnosis: depth_too_shallow_degenerate_per_instance_tests\n");
    else
        printf("diagnosis: no_superlinear_rise; cull time ~linear in N, ns/inst stable\n");

    printf("\n===== shallow H=0 degeneration check =====\n");
    for (int ti = 0; ti < 3; ti++) {
        double step = ti ? rows[ti].shallow_us/rows[ti-1].shallow_us : 1.0;
        double tax = rows[ti].shallow_us / rows[ti].cull_us;
        printf("tier=%-3s n=%4d depth_avg=%.2f inst_tests/N=%.3f "
               "shallow_us=%9.3f step=%5.2f vs_chosen=%.2fx\n",
               LABEL[ti], rows[ti].n, rows[ti].shallow_depth,
               rows[ti].shallow_inst_frac, rows[ti].shallow_us, step, tax);
    }

    /* 三组十二个数: 每组依次 submissions, cull_us, missed, depth_avg */
    printf("\n===== twelve numbers (3 tiers x [submissions cull_us missed depth_avg]) =====\n");
    for (int ti = 0; ti < 3; ti++)
        printf("ROW_%s n=%d submissions=%ld cull_us=%.3f missed=%ld depth_avg=%.2f\n",
               LABEL[ti], rows[ti].n, rows[ti].subs, rows[ti].cull_us,
               rows[ti].miss, rows[ti].depth_avg);
    printf("TWELVE: %ld,%.3f,%ld,%.2f; %ld,%.3f,%ld,%.2f; %ld,%.3f,%ld,%.2f\n",
           rows[0].subs, rows[0].cull_us, rows[0].miss, rows[0].depth_avg,
           rows[1].subs, rows[1].cull_us, rows[1].miss, rows[1].depth_avg,
           rows[2].subs, rows[2].cull_us, rows[2].miss, rows[2].depth_avg);
    printf("METRICS_JSON: {\"tiers\":[200,400,800],\"k_ns_per_submit\":%.1f,"
           "\"submissions\":[%ld,%ld,%ld],\"cull_us\":[%.3f,%.3f,%.3f],"
           "\"missed\":[%ld,%ld,%ld],\"depth_avg\":[%.2f,%.2f,%.2f],"
           "\"gran_avg\":[%.2f,%.2f,%.2f],\"hstar\":[%d,%d,%d]}\n",
           K_NS_PER_SUBMIT,
           rows[0].subs, rows[1].subs, rows[2].subs,
           rows[0].cull_us, rows[1].cull_us, rows[2].cull_us,
           rows[0].miss, rows[1].miss, rows[2].miss,
           rows[0].depth_avg, rows[1].depth_avg, rows[2].depth_avg,
           rows[0].gran_avg, rows[1].gran_avg, rows[2].gran_avg,
           rows[0].hstar, rows[1].hstar, rows[2].hstar);

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
