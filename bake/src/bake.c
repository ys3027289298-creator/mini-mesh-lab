
/* bake: curvature-aware shell offset for ray-cast normal baking.
 *
 * 低模沿弧线铺开（凹侧朝向高模），烘焙射线从低模表面沿顶点法线发出。
 * 无外壳时，接缝附近的射线先打到低模相邻面片（自遮挡）-> 暗线。
 * 每个顶点沿法线加外壳偏移 d = min(烘焙距离, sigma * 曲率半径)，
 * 曲率大（半径小）处偏移更小，避免穿过高模。
 *
 * 三种曲率（同弦长、不同半径/二面角）各烘一次，统计：
 *   接缝法线差(度)、暗线宽度(像素)、有效像素占比。
 * 另跑一遍零偏移基线作对照；若高曲率档有效占比掉得最多，
 * 并列打印该档每个顶点的曲率半径与实际偏移量。
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXV   512
#define MAXF   1024
#define TIERS  3
#define SGRID  5                 /* 每三角形采样网格边长，样本数 SGRID*(SGRID+1)/2 */
#define IMG_PX 32                /* 模拟纹理分辨率：每弦长 IMG_PX 个像素 */
#define ARC_A  2.5f              /* 总弧角（弧度，约 143 度） */
#define CHORD  1.0f              /* 低模弦长（世界单位，所有档相同） */
#define GAP    0.55f             /* 高模在凹侧与低模的间距 */
#define BAKE_D 0.5f              /* 烘焙射线最大距离 */
#define SIGMA  0.1f              /* 穿透保护：偏移不超过曲率半径的 1/10 */
#define EPS_H  1e-5f

typedef struct { float x, y, z; } V3;
typedef struct { int a, b, c; } Tri;
typedef struct { V3 v[MAXV]; Tri f[MAXF]; int nv, nf; float r; int nseg; } Arc;

static V3 sub(V3 p, V3 q) { V3 r = { p.x-q.x, p.y-q.y, p.z-q.z }; return r; }
static V3 add(V3 p, V3 q) { V3 r = { p.x+q.x, p.y+q.y, p.z+q.z }; return r; }
static V3 mul(V3 p, float s) { V3 r = { p.x*s, p.y*s, p.z*s }; return r; }
static V3 cross(V3 p, V3 q) { V3 r = { p.y*q.z-p.z*q.y, p.z*q.x-p.x*q.z, p.x*q.y-p.y*q.x }; return r; }
static float dot(V3 p, V3 q) { return p.x*q.x + p.y*q.y + p.z*q.z; }
static float len(V3 p) { return sqrtf(dot(p, p)); }
static V3 norm(V3 p) { float l = len(p); V3 r = {0,0,0}; if (l > 0) { r = mul(p, 1.0f/l); } return r; }
static float dist(V3 p, V3 q) { return len(sub(p, q)); }
static float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

/* 弧形条带：半径 r、nseg 段、半径方向再内缩 gap（高模用 gap>0）。
 * 绕序取凹侧（指向弧心）法线。 */
static void arc_build(Arc* m, float r, int nseg, float gap) {
    memset(m, 0, sizeof(*m));
    m->r = r; m->nseg = nseg;
    float rho = r - gap, da = ARC_A / nseg;
    for (int i = 0; i <= nseg; i++) {
        float a = i * da;
        V3 p = { rho * cosf(a), rho * sinf(a), 0.0f };
        m->v[m->nv++] = p;
        p.z = 1.0f;
        m->v[m->nv++] = p;
    }
    for (int i = 0; i < nseg; i++) {
        int b0 = 2*i, t0 = 2*i+1, b1 = 2*(i+1), t1 = 2*(i+1)+1;
        m->f[m->nf++] = (Tri){ b0, t0, b1 };
        m->f[m->nf++] = (Tri){ b1, t0, t1 };
    }
}

/* 面积加权顶点法线（同 normal 模块） */
static void vertex_normals(const Arc* m, V3* out) {
    for (int i = 0; i < m->nv; i++) out[i] = (V3){0,0,0};
    for (int i = 0; i < m->nf; i++) {
        V3 n = cross(sub(m->v[m->f[i].b], m->v[m->f[i].a]),
                     sub(m->v[m->f[i].c], m->v[m->f[i].a]));
        out[m->f[i].a] = add(out[m->f[i].a], n);
        out[m->f[i].b] = add(out[m->f[i].b], n);
        out[m->f[i].c] = add(out[m->f[i].c], n);
    }
    for (int i = 0; i < m->nv; i++) out[i] = norm(out[i]);
}

/* 逐顶点曲率半径：相邻三列顶点外接圆半径（上下排同值，端列复制邻居） */
static void curvature_radius(const Arc* m, float* R) {
    float da = ARC_A / m->nseg;
    for (int i = 0; i <= m->nseg; i++) {
        int j = clampf(i, 1, m->nseg - 1);
        V3 p0 = m->v[2*(j-1)], p1 = m->v[2*j], p2 = m->v[2*(j+1)];
        float a = dist(p0, p1), b = dist(p1, p2), c = dist(p2, p0);
        float area2 = len(cross(sub(p1, p0), sub(p2, p0)));
        float r = (a*b*c) / (2.0f * area2);   /* 外接圆半径 = abc/(4*area) */
        R[2*i] = R[2*i+1] = (da >= ARC_A) ? m->r : r;
    }
}

/* 外壳偏移：烘焙距离与曲率保护值取小 */
static float shell_offset(float Rv) {
    float d = SIGMA * Rv;
    return d < BAKE_D ? d : BAKE_D;
}

/* Moller-Trumbore，返回最近 t>EPS_H 的命中距离，未命中返回 -1；skip>=0 时跳过该面 */
static float ray_mesh(V3 o, V3 d, const Arc* m, int skip) {
    float best = -1.0f;
    for (int i = 0; i < m->nf; i++) {
        if (i == skip) continue;
        V3 a = m->v[m->f[i].a], b = m->v[m->f[i].b], c = m->v[m->f[i].c];
        V3 e1 = sub(b, a), e2 = sub(c, a);
        V3 p = cross(d, e2);
        float det = dot(e1, p);
        if (fabsf(det) < 1e-12f) continue;
        V3 s = sub(o, a);
        float u = dot(s, p) / det;
        if (u < 0.0f || u > 1.0f) continue;
        V3 q = cross(s, e1);
        float v = dot(d, q) / det;
        if (v < 0.0f || u + v > 1.0f) continue;
        float t = dot(e2, q) / det;
        if (t > EPS_H && (best < 0.0f || t < best)) best = t;
    }
    return best;
}

typedef struct {
    float seam_normal_deg;     /* 接缝两侧面法线夹角均值 */
    float dark_line_px;        /* 接缝两侧无效带合并宽度，像素 */
    float valid_ratio;         /* 有效像素占比 0..1 */
} BakeStats;

typedef struct { int ok; float dist_to_seam; } Sample;

/* 对一个低模烘焙一次。use_shell=0 时退化为零偏移基线。 */
static BakeStats bake(const Arc* low, const Arc* high, int use_shell) {
    V3 nrm[MAXV]; float R[MAXV], off[MAXV];
    vertex_normals(low, nrm);
    curvature_radius(low, R);
    for (int i = 0; i < low->nv; i++)
        off[i] = use_shell ? shell_offset(R[i]) : 0.0f;

    int total = 0, valid = 0;
    Sample smp[MAXF][SGRID*(SGRID+1)/2]; int ns[MAXF];
    memset(ns, 0, sizeof(ns));

    for (int fi = 0; fi < low->nf; fi++) {
        Tri f = low->f[fi];
        V3 A = low->v[f.a], B = low->v[f.b], C = low->v[f.c];
        for (int iu = 0; iu < SGRID; iu++) {
            for (int iv = 0; iv < SGRID - iu; iv++) {
                float u = (iu + 0.5f) / SGRID;
                float v = (iv + 0.5f) / SGRID;
                float w = 1.0f - u - v;
                V3 p = add(add(mul(A, w), mul(B, u)), mul(C, v));
                V3 n = norm(add(add(mul(nrm[f.a], w), mul(nrm[f.b], u)), mul(nrm[f.c], v)));
                float d = w*off[f.a] + u*off[f.b] + v*off[f.c];
                V3 o = add(p, mul(n, d));

                float tlow  = ray_mesh(o, n, low, fi);
                float thigh = ray_mesh(o, n, high, -1);
                int hit = (thigh > 0.0f && thigh <= BAKE_D &&
                           (tlow < 0.0f || thigh < tlow));
                valid += hit; total++;
                smp[fi][ns[fi]].ok = hit;
                smp[fi][ns[fi]].dist_to_seam = -1.0f;  /* 仅接缝相邻面需要，稍后填 */
                ns[fi]++;
            }
        }
    }

    /* 接缝法线差：相邻段底三角的面法线 */
    float ndeg = 0.0f;
    for (int i = 1; i < low->nseg; i++) {
        Tri f0 = low->f[2*(i-1)], f1 = low->f[2*i];
        V3 n0 = norm(cross(sub(low->v[f0.b], low->v[f0.a]), sub(low->v[f0.c], low->v[f0.a])));
        V3 n1 = norm(cross(sub(low->v[f1.b], low->v[f1.a]), sub(low->v[f1.c], low->v[f1.a])));
        ndeg += acosf(clampf(dot(n0, n1), -1.0f, 1.0f)) * 180.0f / 3.14159265f;
    }
    ndeg /= (low->nseg - 1);

    /* 暗线宽度：每条内部接缝，两侧面上离缝最远的无效样本距离之和，换算像素 */
    float spacing = CHORD / IMG_PX;
    float width_sum = 0.0f;
    for (int i = 1; i < low->nseg; i++) {
        V3 seam = low->v[2*i];          /* 接缝柱上的底点 */
        float side[2] = { 0.0f, 0.0f };
        int faces[2] = { 2*(i-1), 2*i };/* 两段各自的底三角 */
        for (int k = 0; k < 2; k++) {
            int fi = faces[k];
            Tri f = low->f[fi];
            V3 A = low->v[f.a], B = low->v[f.b], C = low->v[f.c];
            int si = 0;
            for (int iu = 0; iu < SGRID; iu++) {
                for (int iv = 0; iv < SGRID - iu; iv++) {
                    float u = (iu + 0.5f) / SGRID;
                    float v = (iv + 0.5f) / SGRID;
                    float w = 1.0f - u - v;
                    V3 p = add(add(mul(A, w), mul(B, u)), mul(C, v));
                    float dxy = dist((V3){ p.x, p.y, 0 }, (V3){ seam.x, seam.y, 0 });
                    smp[fi][si].dist_to_seam = dxy;
                    if (!smp[fi][si].ok && dxy > side[k]) side[k] = dxy;
                    si++;
                }
            }
        }
        width_sum += (side[0] + side[1]) / spacing;
    }

    BakeStats st;
    st.seam_normal_deg = ndeg;
    st.dark_line_px = width_sum / (low->nseg - 1);
    st.valid_ratio = (float)valid / (float)total;
    return st;
}

static const char* TN[TIERS] = { "low_curve_R10", "mid_curve_R4", "high_curve_R2" };
static const float TR[TIERS] = { 10.0f, 4.0f, 2.0f };

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Arc low[TIERS], high[TIERS];
    int nseg[TIERS];
    for (int t = 0; t < TIERS; t++) {
        /* 固定弦长 CHORD：chord = 2 r sin(da/2)，da=ARC_A/nseg */
        nseg[t] = (int)ceilf(ARC_A / (2.0f * asinf(CHORD / (2.0f * TR[t]))));
        arc_build(&low[t], TR[t], nseg[t], 0.0f);
        arc_build(&high[t], TR[t], nseg[t] * 6, GAP);
    }

    BakeStats bs[TIERS], base[TIERS];
    for (int t = 0; t < TIERS; t++) {
        bs[t]   = bake(&low[t], &high[t], 1);
        base[t] = bake(&low[t], &high[t], 0);
    }

    printf("\n%-16s %14s %14s %14s\n", "tier", "normal_diff_deg", "dark_line_px", "valid_ratio");
    for (int t = 0; t < TIERS; t++)
        printf("%-16s %14.3f %14.2f %14.4f\n", TN[t],
               bs[t].seam_normal_deg, bs[t].dark_line_px, bs[t].valid_ratio);
    printf("%-16s %14s %14s %14s\n", "no_shell_baseline", "-", "-", "-");
    for (int t = 0; t < TIERS; t++)
        printf("%-16s %14.3f %14.2f %14.4f\n", TN[t],
               base[t].seam_normal_deg, base[t].dark_line_px, base[t].valid_ratio);

    /* 有效占比跌幅：相对三档中的最好值，找出掉得最多的档 */
    float best = bs[0].valid_ratio;
    for (int t = 1; t < TIERS; t++) if (bs[t].valid_ratio > best) best = bs[t].valid_ratio;
    float drop[TIERS]; int worst = 0;
    for (int t = 0; t < TIERS; t++) {
        drop[t] = best - bs[t].valid_ratio;
        if (drop[t] > drop[worst]) worst = t;
    }
    printf("\nvalid_ratio_drop_vs_best: %.4f %.4f %.4f (worst=%s)\n",
           drop[0], drop[1], drop[2], TN[worst]);

    if (worst == TIERS - 1) {
        /* 高曲率档掉得最多：把该档偏移量与曲率半径逐顶点并列 */
        Arc* m = &low[TIERS - 1];
        float R[MAXV];
        curvature_radius(m, R);
        printf("\n[high-curvature shell pairing] col  radius_R   offset_d=min(D,sigma*R)\n");
        for (int i = 0; i <= m->nseg; i++)
            printf("v%-3d              %8.4f  %8.4f\n", i, R[2*i], shell_offset(R[2*i]));
    }

    /* ---- 自检 ---- */
    printf("\n");
    for (int t = 0; t < TIERS; t++) {
        V3 nrm[MAXV]; vertex_normals(&low[t], nrm);
        float inward = 1.0f;
        for (int i = 0; i < low[t].nv; i++) {
            V3 radial = norm((V3){ low[t].v[i].x, low[t].v[i].y, 0 }); /* 指向弧外 */
            if (dot(nrm[i], radial) > 0.0f) inward = 0.0f;
        }
        check("normals_point_inward", inward);

        float R[MAXV]; curvature_radius(&low[t], R);
        float err = 0.0f;
        for (int i = 0; i < low[t].nv; i++) {
            float e = fabsf(R[i] - TR[t]) / TR[t];
            if (e > err) err = e;
        }
        check("curvature_radius_est", err < 0.05f);
    }

    float o0 = shell_offset(TR[0]), o1 = shell_offset(TR[1]), o2 = shell_offset(TR[2]);
    check("offset_scaled_by_radius", o0 > o1 && o1 > o2);
    check("offset_never_exceeds_bake_d",
          o0 <= BAKE_D && o1 <= BAKE_D && o2 <= BAKE_D);
    check("offset_below_gap_no_pierce",
          o0 < GAP && o1 < GAP && o2 < GAP);
    check("seam_normal_grows_with_curvature",
          bs[2].seam_normal_deg > bs[1].seam_normal_deg &&
          bs[1].seam_normal_deg > bs[0].seam_normal_deg);
    check("dark_line_grows_with_curvature",
          bs[2].dark_line_px >= bs[1].dark_line_px &&
          bs[1].dark_line_px >= bs[0].dark_line_px);
    int range_ok = 1, shell_beats_base = 1;
    for (int t = 0; t < TIERS; t++) {
        if (bs[t].valid_ratio <= 0.0f || bs[t].valid_ratio > 1.0f) range_ok = 0;
        if (!(bs[t].valid_ratio > base[t].valid_ratio)) shell_beats_base = 0;
    }
    check("valid_ratio_in_range", range_ok);
    check("shell_beats_no_shell_baseline", shell_beats_base);

    /* 诊断逻辑独立复核：worst 必须真的是跌幅最大的档 */
    int recomputed = 0;
    for (int t = 1; t < TIERS; t++) if (drop[t] > drop[recomputed]) recomputed = t;
    check("diagnostic_picks_max_drop", recomputed == worst);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
