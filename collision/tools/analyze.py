#!/usr/bin/env python3
"""Analyze collision-pipeline stress data: fit scaling exponents, locate the
knee where proxy (linear) hands over to broadphase (quadratic), and find where
the 16.7 ms frame budget is crossed. Run with: python3 -I analyze.py <csv...>"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["Noto Sans CJK SC", "WenQuanYi Zen Hei", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

# validated categorical palette (dataviz skill, light mode)
BLUE, ORANGE, AQUA, INK = "#2a78d6", "#eb6834", "#1baf7a", "#52514e"
SURFACE, GRID = "#fcfcfb", "#e3e2de"
BUDGET_MS = 1000.0 / 60.0   # 60 FPS frame budget

def load(path):
    rows = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows[int(r["n"])] = {k: float(r[k]) for k in
                                 ("total_ns", "proxy_ns", "broad_ns",
                                  "narrow_ns", "motion_ns", "pairs")}
    return rows

def main(csv_paths, out_dir):
    runs = [load(p) for p in csv_paths]
    ns = np.array(sorted(runs[0]), dtype=float)
    def med(key):
        return np.array([np.median([r[int(n)][key] for r in runs]) for n in ns])

    total, proxy, broad = med("total_ns"), med("proxy_ns"), med("broad_ns")
    narrow, motion = med("narrow_ns"), med("motion_ns")
    total_ms, proxy_ms, broad_ms = total/1e6, proxy/1e6, broad/1e6

    ln = np.log(ns)
    def local_exponent(y):
        ly = np.log(np.maximum(y, 1e-12))
        return np.gradient(ly, ln)
    e_total, e_proxy, e_broad = local_exponent(total), local_exponent(proxy), local_exponent(broad)

    # global power-law fits per phase
    def fit_exp(x, y, lo=None, hi=None):
        m = np.ones(len(x), bool)
        if lo is not None: m &= x >= lo
        if hi is not None: m &= x <= hi
        slope, intercept = np.polyfit(ln[m], np.log(y[m]), 1)
        return slope, np.exp(intercept)
    e_proxy_fit, _ = fit_exp(ns, proxy, None, 220)
    e_broad_fit, _ = fit_exp(ns, broad, 460, None)
    e_broad_asym, _ = fit_exp(ns, broad, 900, None)   # asymptotic window
    e_total_pre, _  = fit_exp(ns, total, None, 320)
    e_total_post, _ = fit_exp(ns, total, 460, None)

    # knee of the total curve, two independent detectors:
    # (a) first n where the noisy local exponent crosses 1.5
    knee_idx = int(np.argmax(e_total >= 1.5))
    knee_n_cross = int(ns[knee_idx])
    # (b) best two-segment piecewise-linear fit on log-log (minimum SSE breakpoint)
    ly_total = np.log(total)
    def piecewise_sse(break_i):
        sse = 0.0
        for lo, hi in ((0, break_i + 1), (break_i, len(ns))):
            if hi - lo < 2: return np.inf
            p = np.polyfit(ln[lo:hi], ly_total[lo:hi], 1)
            sse += float(np.sum((ly_total[lo:hi] - np.polyval(p, ln[lo:hi])) ** 2))
        return sse
    bp = min(range(2, len(ns) - 2), key=piecewise_sse)
    knee_n_fit = int(ns[bp])
    knee_n = knee_n_fit

    # proxy/broad dominance crossover, interpolated between bracketing samples
    diff = proxy - broad
    sign_change = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]
    if len(sign_change):
        i = sign_change[0]
        cross_n = int(round(ns[i] - diff[i] * (ns[i+1] - ns[i]) / (diff[i+1] - diff[i])))
    else:
        cross_n = int(ns[int(np.argmax(broad >= proxy))])

    # frame-budget crossing: median over every bench_ext*.csv next to the inputs
    budget_n_meas = None
    ext_files = sorted(Path(out_dir).glob("bench_ext*.csv"))
    ext_rows = {}
    ext_ns, ext_tot = np.array([]), np.array([])
    for ep in ext_files:
        with open(ep, newline="") as f:
            for r in csv.DictReader(f):
                ext_rows.setdefault(int(r["n"]), []).append(float(r["total_ns"]) / 1e6)
    if ext_rows:
        ext_ns = np.array(sorted(ext_rows), dtype=float)
        ext_tot = np.array([np.median(ext_rows[int(n)]) for n in ext_ns])
        above = np.where(ext_tot > BUDGET_MS)[0]
        if len(above):
            i = above[0]
            if i > 0:
                x0, x1, y0, y1 = ext_ns[i-1], ext_ns[i], ext_tot[i-1], ext_tot[i]
                budget_n_meas = x0 + (BUDGET_MS-y0)*(x1-x0)/(y1-y0)

    # analytic fallback/estimate: fit total = a*n + b*n^2 over the full range
    A = np.column_stack([ns, ns**2])
    coef, *_ = np.linalg.lstsq(A, total_ms, rcond=None)
    a, b = coef
    roots = np.roots([b, a, -BUDGET_MS])
    budget_n_fit = float(max(roots))

    # cache-boundary accounting (machine: L1d 48 KiB, L2 1 MiB, L3 32 MiB)
    hdr = Path(__file__).resolve().parent.parent / "src" / "collision.h"
    bv = 256
    for line in hdr.read_text().splitlines():
        if line.startswith("#define BODY_VERTS"):
            bv = int(line.split()[2])
    vert_bytes = bv * 3 * 4          # rotated point cloud per object
    proxy_bytes = 6 * 4              # one AABB = 6 floats
    caches = [("L1d", 48*1024), ("L2", 1024*1024), ("L3", 32*1024*1024)]
    cache_info = {}
    for name, cap in caches:
        cache_info[name] = {
            "capacity_KiB": cap // 1024,
            "verts_pool_fits_until_n": cap // vert_bytes,
            "proxy_array_fits_until_n": cap // proxy_bytes,
        }

    # phase share table
    share = [(int(n), 100*p/t, 100*br/t) for n, t, p, br in
             zip(ns, total_ms, proxy_ms, broad_ms)]

    results = {
        "runs": [str(p) for p in csv_paths],
        "exponent_proxy_fit_n_le_220": round(e_proxy_fit, 3),
        "exponent_broad_fit_n_ge_460": round(e_broad_fit, 3),
        "exponent_broad_asymptotic_n_ge_900": round(e_broad_asym, 3),
        "exponent_total_pre_knee": round(e_total_pre, 3),
        "exponent_total_post_knee": round(e_total_post, 3),
        "knee_n_piecewise_fit": knee_n_fit,
        "knee_n_local_exp_1p5": knee_n_cross,
        "proxy_broad_crossover_n_interp": cross_n,
        "frame_budget_ms": round(BUDGET_MS, 3),
        "max_objects_measured_interp": None if budget_n_meas is None else round(budget_n_meas),
        "max_objects_quadratic_fit": round(budget_n_fit),
        "fit_total_ms_eq": f"{a:.3e}*n + {b:.3e}*n^2",
        "cache_boundaries": cache_info,
        "phase_share_pct": [{"n": n, "proxy": round(p, 1), "broad": round(br, 1)}
                            for n, p, br in share],
        "local_exponents": [{"n": int(n), "e_total": round(float(t), 3),
                             "e_proxy": round(float(p), 3),
                             "e_broad": round(float(bb), 3)}
                            for n, t, p, bb in zip(ns, e_total, e_proxy, e_broad)],
    }
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "analysis.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))

    # ---------- chart 1: log-log phase breakdown ----------
    fig, ax = plt.subplots(figsize=(10, 6.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    # transition band: from first measurable quadratic term to the phase handover
    ax.axvspan(knee_n_fit, cross_n, color="#d9d8d2", alpha=0.45, zorder=0)
    ax.text((knee_n_fit*cross_n)**0.5, 0.0035, "过渡区\n线性主导 → 二次主导",
            ha="center", va="bottom", fontsize=9.5, color=INK)
    ax.plot(ns, total_ms, "-o", color=INK, lw=2, ms=6, label="单帧总耗时", zorder=3)
    ax.plot(ns, proxy_ms, "-s", color=BLUE, lw=2, ms=6, label="代理生成（O(n)）", zorder=3)
    ax.plot(ns, broad_ms, "-^", color=ORANGE, lw=2, ms=6,
            label="粗检测 AABB 两两测试（O(n²)）", zorder=3)
    ax.plot(ns, narrow/1e6, ":d", color=AQUA, lw=2, ms=5,
            label="精检测（候选对，近似 O(n)）", zorder=3)
    if len(ext_ns):
        ax.plot(ext_ns, ext_tot, "o", color=INK, ms=8, mfc="white", mew=2,
                label="外推实测点（>2000）", zorder=4)
    # reference slopes: same color family, direct annotated (no legend entry)
    ref_n = np.array([40, 200]); ax.plot(ref_n, proxy_ms[3]*ref_n/ns[3], "--",
                                         color=BLUE, lw=1.2, alpha=0.5, zorder=2)
    ax.text(46, proxy_ms[3]*46/ns[3]*0.72, "斜率 1", color=BLUE, fontsize=9, alpha=0.85)
    ref_n2 = np.array([460, 2000]); ax.plot(ref_n2, broad_ms[8]*(ref_n2/ns[8])**2,
                                            "--", color=ORANGE, lw=1.2, alpha=0.5, zorder=2)
    ax.text(1050, broad_ms[8]*(1050/ns[8])**2*1.18, "斜率 2", color=ORANGE,
            fontsize=9, alpha=0.85)
    ax.axhline(BUDGET_MS, color="#e34948", lw=1.6, ls="-.", zorder=2)
    ax.text(12, BUDGET_MS*1.08, "16.7 ms 帧预算（60 FPS）", color="#e34948", fontsize=10)
    ax.axvline(cross_n, color="#0f7d55", lw=1.4, ls=":", alpha=0.9, zorder=2)
    ax.text(cross_n*1.03, 0.35, f"瓶颈易主 n≈{cross_n}\n（粗检测耗时追平代理生成）",
            color="#0f7d55", fontsize=9.5)
    ax.set_xscale("log", base=10); ax.set_yscale("log", base=10)
    ax.set_xlabel("物体数 n", fontsize=11)
    ax.set_ylabel("单帧耗时 (ms, 对数轴)", fontsize=11)
    ax.set_title("碰撞检测流水线压测：单帧耗时随物体数的增长（10 → 2000，密度固定）",
                 fontsize=12, pad=12)
    ax.grid(True, which="both", color=GRID, lw=0.7)
    ax.legend(fontsize=9.5, frameon=False, loc="lower right")
    for sp in ax.spines.values(): sp.set_color("#c9c8c2")
    fig.tight_layout()
    fig.savefig(out / "curve_loglog.png", dpi=150)

    # ---------- chart 2: local exponents (3-point smoothed) ----------
    def smooth(x):
        x = np.asarray(x, float)
        out = x.copy()
        out[1:-1] = (x[:-2] + 2*x[1:-1] + x[2:]) / 4.0
        return out
    s_total, s_proxy, s_broad = smooth(e_total), smooth(e_proxy), smooth(e_broad)
    fig2, ax2 = plt.subplots(figsize=(10, 5.4), facecolor=SURFACE)
    ax2.set_facecolor(SURFACE)
    ax2.axvspan(knee_n_fit, cross_n, color="#d9d8d2", alpha=0.45, zorder=0)
    ax2.plot(ns, s_total, "-o", color=INK, lw=2, ms=6, label="总耗时局部阶数（平滑）", zorder=3)
    ax2.plot(ns, s_proxy, "-s", color=BLUE, lw=2, ms=6, label="代理生成局部阶数（平滑）", zorder=3)
    ax2.plot(ns, s_broad, "-^", color=ORANGE, lw=2, ms=6,
             label="粗检测局部阶数（平滑，小 n 受计时噪声影响）", zorder=3)
    ax2.axhline(1.0, color=BLUE, lw=1, ls="--", alpha=0.5)
    ax2.axhline(2.0, color=ORANGE, lw=1, ls="--", alpha=0.5)
    ax2.axhline(1.5, color="#e34948", lw=1.2, ls="-.")
    ax2.text(11, 1.56, "阶数 1.5（拐点判据）", color="#e34948", fontsize=10)
    ax2.axvline(cross_n, color="#0f7d55", lw=1.4, ls=":", alpha=0.9)
    ax2.text(cross_n*1.03, 0.35, f"总阶数越过 1.5\nn≈{knee_n_cross}", color="#0f7d55", fontsize=9.5)
    ax2.set_xscale("log", base=10)
    ax2.set_ylim(0, 2.6)
    ax2.set_xlabel("物体数 n", fontsize=11)
    ax2.set_ylabel("局部增长阶数 d(log t) / d(log n)", fontsize=11)
    ax2.set_title("局部阶数曲线：线性段（≈1）向二次段（≈2）的切换", fontsize=12, pad=12)
    ax2.grid(True, which="both", color=GRID, lw=0.7)
    ax2.legend(fontsize=9.5, frameon=False, loc="upper left")
    for sp in ax2.spines.values(): sp.set_color("#c9c8c2")
    fig2.tight_layout()
    fig2.savefig(out / "exponents.png", dpi=150)

    print(json.dumps(results, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    args = sys.argv[1:]
    csvs = [a for a in args if not a.startswith("--")]
    out_dir = str(Path(csvs[0]).resolve().parent)
    main(csvs, out_dir)
