#!/usr/bin/env python3
"""Vertex-color quantization: banding vs. storage growth, bit-depth inversion.

Pipeline per tier:
  1. build a smooth neutral-axis sRGB gradient spanning a local color range;
  2. quantize each channel to b bits (round-to-nearest, clamp), dequantize;
  3. measure bands (distinct steps), max single-step dE00, round-trip error;
  4. storage growth under two storage models.

All perceptual color differences use CIEDE2000 on sRGB-derived CIELAB
(D65). Run with -I; depends only on numpy + matplotlib.
"""
import csv
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, gridspec

# ---------------------------------------------------------------- constants
# sRGB D65 white point (CIE 1931, normalized)
D65 = np.array([0.95047, 1.0, 1.08883])

# Three local color-difference ranges. Defined as per-channel sRGB span on a
# neutral (r=g=b) ramp; the delta-E span is measured, not assumed.
# Subtle: 8/255 (~3.1% channel span) — e.g. baked AO / soft gradient.
# Medium: 64/255 (~25%).
# Full:   255/255 (the whole 0..1 gamut) — bit depth chosen against this.
TIERS = [
    ("narrow 窄范围", 8 / 255.0),
    ("medium 中范围", 64 / 255.0),
    ("full 全色域", 255 / 255.0),
]

BITS_MIN, BITS_MAX = 8, 16
# Acceptable single-step for a SMOOTH GRADIENT. The classic ~1 dE00 is a
# just-noticeable difference for two isolated patches, but a Mach-band edge
# in a continuous gradient makes steps roughly twice as visible; the display
# engineering no-contouring criterion is therefore ~0.5 JND = 0.2 dE00.
STEP_TOL = 0.2
N_SAMPLES = 70001  # dense enough to resolve all 65535 codes of 16-bit full range

# dataviz reference palette (validated)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, SURFACE = "#e1e0d9", "#fcfcfb"

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output")


# ------------------------------------------------------------- color science
def srgb_to_lin(c):
    c = np.asarray(c, dtype=np.float64)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def lin_to_xyz(rgb_lin):
    m = np.array(
        [
            [0.4123908, 0.3575843, 0.1804808],
            [0.2126390, 0.7151687, 0.0721923],
            [0.0193308, 0.1191948, 0.9505322],
        ]
    )
    return rgb_lin @ m.T


def xyz_to_lab(xyz):
    xyz = np.asarray(xyz, dtype=np.float64) / D65
    eps = 216 / 24389.0
    kap = 24389 / 27.0
    f = np.where(xyz > eps, np.cbrt(xyz), (kap * xyz + 16) / 116)
    L = 116 * f[..., 1] - 16
    a = 500 * (f[..., 0] - f[..., 1])
    b = 200 * (f[..., 1] - f[..., 2])
    return np.stack([L, a, b], axis=-1)


def srgb_to_lab(rgb):
    return xyz_to_lab(lin_to_xyz(srgb_to_lin(rgb)))


def ciede2000(lab1, lab2):
    """CIEDE2000 between two (..., 3) Lab arrays. Sharma et al. 2005."""
    lab1 = np.asarray(lab1, dtype=np.float64)
    lab2 = np.asarray(lab2, dtype=np.float64)
    L1, a1, b1 = lab1[..., 0], lab1[..., 1], lab1[..., 2]
    L2, a2, b2 = lab2[..., 0], lab2[..., 1], lab2[..., 2]

    C1 = np.hypot(a1, b1)
    C2 = np.hypot(a2, b2)
    Cbar = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cbar**7 / (Cbar**7 + 25**7)))
    ap1, ap2 = (1 + G) * a1, (1 + G) * a2
    Cp1, Cp2 = np.hypot(ap1, b1), np.hypot(ap2, b2)
    # arctan2 gives (-180, 180]; mod 360 maps negative angles into [0, 360)
    hp1 = np.degrees(np.arctan2(b1, ap1)) % 360.0
    hp2 = np.degrees(np.arctan2(b2, ap2)) % 360.0

    dLp = L2 - L1
    dCp = Cp2 - Cp1
    dhp = np.where(
        np.abs(hp1 - hp2) <= 180,
        hp2 - hp1,
        hp2 - hp1 + np.where(hp2 <= hp1, 360, -360),
    )
    dhp = np.where(Cp1 * Cp2 == 0, 0.0, dhp)
    dHp = 2 * np.sqrt(Cp1 * Cp2) * np.sin(np.radians(dhp / 2))

    Lbp = (L1 + L2) / 2
    Cbp = (Cp1 + Cp2) / 2
    hbp = np.where(
        Cp1 * Cp2 == 0,
        hp1 + hp2,
        np.where(
            np.abs(hp1 - hp2) <= 180,
            (hp1 + hp2) / 2,
            (hp1 + hp2 + np.where(hp1 + hp2 < 360, 360, -360)) / 2,
        ),
    )
    T = (
        1
        - 0.17 * np.cos(np.radians(hbp - 30))
        + 0.24 * np.cos(np.radians(2 * hbp))
        + 0.32 * np.cos(np.radians(3 * hbp + 6))
        - 0.20 * np.cos(np.radians(4 * hbp - 63))
    )
    dtheta = 30 * np.exp(-(((hbp - 275) / 25) ** 2))
    RC = 2 * np.sqrt(Cbp**7 / (Cbp**7 + 25**7))
    SL = 1 + (0.015 * (Lbp - 50) ** 2) / np.sqrt(20 + (Lbp - 50) ** 2)
    SC = 1 + 0.045 * Cbp
    SH = 1 + 0.015 * Cbp * T
    RT = -np.sin(np.radians(2 * dtheta)) * RC

    kL = kC = kH = 1.0
    return np.sqrt(
        (dLp / (kL * SL)) ** 2
        + (dCp / (kC * SC)) ** 2
        + (dHp / (kH * SH)) ** 2
        + RT * (dCp / (kC * SC)) * (dHp / (kH * SH))
    )


# ------------------------------------------------------------- quantization
def quantize_roundtrip(rgb01, bits):
    """Per-channel round-to-nearest quantization then dequantization."""
    levels = 2**bits - 1
    codes = np.round(np.clip(rgb01, 0, 1) * levels)
    return codes / levels, codes


def storage_growth_doubling(bits):
    """User's model: storage type doubles per added bit. 8-bit = 1x."""
    return 2.0 ** (bits - 8)


def storage_growth_packed(bits):
    """Practical packed channels: rounded up to 8/16/32-bit container widths."""
    width = 8 if bits <= 8 else 16 if bits <= 16 else 32
    return width / 8.0


# ------------------------------------------------------------------- ramp
def make_ramp(span):
    """Smooth neutral-axis sRGB ramp whose per-channel span equals `span`."""
    t = np.linspace(0.0, 1.0, N_SAMPLES)[:, None]
    lo = (1.0 - span) / 2.0  # center the range in the gamut
    return np.repeat(lo + span * t, 3, axis=1)


def ramp_metrics(span, bits):
    rgb = make_ramp(span)
    rgb_q, codes = quantize_roundtrip(rgb, bits)

    lab0, labq = srgb_to_lab(rgb), srgb_to_lab(rgb_q)
    # bands: distinct quantized levels along the ramp
    unique_codes = np.unique(codes[:, 0])
    n_bands = int(len(unique_codes) - 1)
    # max single-step color difference between adjacent surviving levels
    level_rgb = (unique_codes / (2**bits - 1))[:, None].repeat(3, axis=1)
    level_lab = srgb_to_lab(level_rgb)
    step_de = ciede2000(level_lab[:-1], level_lab[1:])
    max_step_de = float(step_de.max()) if step_de.size else 0.0
    # round-trip color error over the whole sampled ramp
    rt_de = ciede2000(lab0, labq)
    return {
        "bands": n_bands,
        "max_step_de": max_step_de,
        "rt_max_de": float(rt_de.max()),
        "rt_rms_de": float(np.sqrt((rt_de**2).mean())),
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    bits = np.arange(BITS_MIN, BITS_MAX + 1)

    # ---- curve 1 & 2 over every bit depth, per tier
    rows = []
    curves = {name: {"bands": [], "step": [], "rtmax": [], "rtrms": []}
              for name, _ in TIERS}
    for name, span in TIERS:
        for b in bits:
            m = ramp_metrics(span, int(b))
            rows.append([name, span, int(b), m["bands"],
                         round(m["max_step_de"], 4),
                         round(m["rt_max_de"], 4), round(m["rt_rms_de"], 4),
                         f"{storage_growth_doubling(int(b)):.1f}",
                         f"{storage_growth_packed(int(b)):.1f}"])
            curves[name]["bands"].append(m["bands"])
            curves[name]["step"].append(m["max_step_de"])
            curves[name]["rtmax"].append(m["rt_max_de"])
            curves[name]["rtrms"].append(m["rt_rms_de"])

    with open(os.path.join(OUT, "quant_curves.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tier", "channel_span", "bits", "bands",
                    "max_step_dE00", "rt_max_dE00", "rt_rms_dE00",
                    "vol_x_doubling_model", "vol_x_packed_model"])
        w.writerows(rows)

    # ---- invert bit depth: two criteria, both reported
    full_step = np.array(curves[TIERS[-1][0]]["step"])
    # criterion 1: per-step perceptual contouring on the full gamut
    chosen = int(next(b for b, d in zip(bits, full_step) if d <= STEP_TOL))
    # classic isolated-patch JND, shown only to explain why it is unsuitable
    b_at_jnd1 = int(next(b for b, d in zip(bits, full_step) if d <= 1.0))
    # criterion 2: every tier keeps at least as many bands as the 8-bit full
    # gamut offers (255) — band count, not step size, governs the narrow tier
    baseline_bands = next(m["bands"] for (n, _), m in
                          zip(TIERS, [ramp_metrics(s, 8) for _, s in TIERS])
                          if n == TIERS[-1][0])
    b_for_bandcount = {}
    for name, span in TIERS:
        b_for_bandcount[name] = int(next(
            b for b, nb in zip(bits, curves[name]["bands"])
            if nb >= baseline_bands))
    chosen_bc = max(b_for_bandcount.values())  # one global depth serving all tiers

    growth_d = storage_growth_doubling(chosen)
    growth_p = storage_growth_packed(chosen)
    print(f"curve 1 inversion (max single-step dE00 <= {STEP_TOL}, "
          f"smooth-gradient half-JND): b = {chosen}")
    print(f"  same criterion with isolated-patch JND=1 gives b = {b_at_jnd1} "
          f"(rejected: 8-bit steps already under 1 yet visibly band)")
    print(f"curve 2 inversion (all tiers >= {baseline_bands} bands, the 8-bit "
          f"full-gamut reference): b = {chosen_bc}")
    for name, _ in TIERS:
        print(f"    {name}: needs {b_for_bandcount[name]} bit locally")
    print(f"  volume at b={chosen}: {growth_d:.1f}x doubling model, "
          f"{growth_p:.1f}x packed; at b={chosen_bc}: "
          f"{storage_growth_doubling(chosen_bc):.1f}x / "
          f"{storage_growth_packed(chosen_bc):.1f}x")

    # ---- the 3 tiers x 4 numbers round-trip at the curve-1 landing depth
    trip = []
    for name, span in TIERS:
        m = ramp_metrics(span, chosen)
        trip.append((name, span, m))
        print(f"{name}: span={span*255:.0f}/255 ch, bands={m['bands']}, "
              f"maxStepDE={m['max_step_de']:.3f}, "
              f"vol x{growth_d:.0f}(doubling)/x{growth_p:.0f}(packed), "
              f"rtMaxDE={m['rt_max_de']:.3f}, rtRmsDE={m['rt_rms_de']:.3f}")

    # local color range vs. depth each tier needs under both criteria
    print("\nlocal range vs. depth (curve-1 step criterion / curve-2 band count):")
    local_need = {}
    for name, span in TIERS:
        steps = np.array(curves[name]["step"])
        bneed = int(next(b for b, d in zip(bits, steps) if d <= STEP_TOL))
        local_need[name] = bneed
        rgb_lo = (1 - span) / 2
        lab_lo = srgb_to_lab(np.full(3, rgb_lo))
        lab_hi = srgb_to_lab(np.full(3, rgb_lo + span))
        span_de = float(ciede2000(lab_lo, lab_hi))
        print(f"  {name}: local span {span*255:3.0f}/255 ch = {span_de:5.1f} dE00, "
              f"step-criterion {bneed} bit, band-count criterion "
              f"{b_for_bandcount[name]} bit, global choice {chosen} bit")

    with open(os.path.join(OUT, "roundtrip_at_chosen.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tier", "channel_span", "bits", "bands", "max_step_dE00",
                    "vol_x_doubling", "vol_x_packed", "rt_max_dE00", "rt_rms_dE00",
                    "bits_for_step_criterion", "bits_for_bandcount_criterion"])
        for (name, span, m) in trip:
            w.writerow([name, span, chosen, m["bands"],
                        round(m["max_step_de"], 4), growth_d, growth_p,
                        round(m["rt_max_de"], 4), round(m["rt_rms_de"], 4),
                        local_need[name], b_for_bandcount[name]])

    make_figure(bits, curves, chosen, chosen_bc, b_at_jnd1, STEP_TOL,
                baseline_bands, b_for_bandcount, local_need, trip)


# ------------------------------------------------------------------- figure
def make_figure(bits, curves, chosen, chosen_bc, b_jnd1, step_tol,
                baseline_bands, b_for_bandcount, local_need, trip):
    for f in ["Noto Sans CJK SC", "WenQuanYi Zen Hei", "DejaVu Sans"]:
        if any(f in fn.name for fn in font_manager.fontManager.ttflist):
            plt.rcParams["font.sans-serif"] = [f]
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["mathtext.fontset"] = "dejavusans"
    plt.rcParams["figure.facecolor"] = SURFACE
    plt.rcParams["axes.facecolor"] = SURFACE
    plt.rcParams["text.color"] = INK
    plt.rcParams["axes.edgecolor"] = "#c3c2b7"
    plt.rcParams["axes.labelcolor"] = INK2
    plt.rcParams["xtick.color"] = MUTED
    plt.rcParams["ytick.color"] = MUTED

    tier_colors = [AQUA, ORANGE, BLUE]
    fig = plt.figure(figsize=(12, 9.5), dpi=150)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.15, 1], hspace=0.36, wspace=0.24)

    # panel A: descending curve — max single-step dE00 vs bits (all tiers)
    axA = fig.add_subplot(gs[0, :])
    for (name, _), col in zip(TIERS, tier_colors):
        axA.plot(bits, curves[name]["step"], marker="o", lw=2, ms=6,
                 color=col, label=name)
    axA.axhline(step_tol, color=MUTED, lw=1.4, ls="--")
    axA.annotate(f"渐变无条带阈值（半 JND）$\\Delta E_{{00}}={step_tol}$",
                 xy=(BITS_MAX - 0.15, step_tol),
                 xytext=(BITS_MAX - 0.15, step_tol * 1.9), ha="right",
                 color=INK2, fontsize=9)
    axA.axhline(1.0, color=MUTED, lw=1.0, ls=(0, (1, 3)))
    axA.annotate("色块 JND = 1（8 位已低于它，不能用作判据）",
                 xy=(BITS_MAX - 0.15, 1.0),
                 xytext=(BITS_MAX - 0.15, 1.45), ha="right",
                 color=MUTED, fontsize=8.5)
    axA.axvline(chosen, color=INK, lw=1.2, ls=":")
    axA.annotate(f"反解落点 b = {chosen}", xy=(chosen, 0.008),
                 xytext=(chosen + 0.12, 0.011), color=INK, fontsize=10)
    axA.set_yscale("log")
    axA.set_xlabel("每通道位数 b (bit)")
    axA.set_ylabel(r"最大单级色差 $\Delta E_{00}$（对数轴）")
    axA.set_title("曲线①  台阶严重度随位数下降（三条为不同局部色差范围）",
                  loc="left", fontsize=12, color=INK, pad=8)
    axA.grid(True, color=GRID, lw=0.8)
    axA.legend(frameon=False, fontsize=9, ncol=3, loc="upper right")

    # panel B: ascending curve — storage growth vs bits, two models
    axB = fig.add_subplot(gs[1, 0])
    vol_d = [storage_growth_doubling(int(b)) for b in bits]
    vol_p = [storage_growth_packed(int(b)) for b in bits]
    axB.plot(bits, vol_d, marker="o", lw=2, ms=6, color=BLUE,
             label="逐位翻倍模型（题设）")
    axB.plot(bits, vol_p, marker="s", lw=2, ms=6, color=ORANGE,
             label="实际打包 8/16/32 位容器")
    axB.set_yscale("log", base=2)
    axB.axvline(chosen, color=INK, lw=1.2, ls=":")
    axB.annotate(f"b={chosen}: ×{storage_growth_doubling(chosen):.0f} / "
                 f"×{storage_growth_packed(chosen):.0f}",
                 xy=(chosen, storage_growth_doubling(chosen)),
                 xytext=(9.4, 55), color=INK, fontsize=9)
    axB.axvline(chosen_bc, color=AQUA, lw=1.2, ls=":")
    axB.annotate(f"台阶数判据 b={chosen_bc}: "
                 f"×{storage_growth_doubling(chosen_bc):.0f} / "
                 f"×{storage_growth_packed(chosen_bc):.0f}",
                 xy=(chosen_bc, 2), xytext=(8.15, 4.5), color="#0d6b4d",
                 fontsize=8.5)
    axB.set_xlim(7.7, 16.3)
    axB.set_xlabel("每通道位数 b (bit)")
    axB.set_ylabel(r"顶点色体积倍数（$\log_2$ 轴）")
    axB.set_title("曲线②  体积增长随位数上升", loc="left", fontsize=12,
                  color=INK, pad=8)
    axB.grid(True, color=GRID, lw=0.8)
    axB.legend(frameon=False, fontsize=9, loc="upper left")

    # panel C: bands at chosen depth per tier — the visible-step landing.
    # Log scale: 17 vs 511 spans 1.5 decades, a linear axis hides the 17 bar.
    axC = fig.add_subplot(gs[1, 1])
    names = [t[0].split()[0] for t in TIERS]
    band_counts = [m["bands"] for _, _, m in trip]
    x = np.arange(len(names))
    bars = axC.bar(x, band_counts, color=tier_colors, width=0.62,
                   edgecolor=SURFACE, linewidth=2)
    axC.set_yscale("log", base=2)
    axC.set_xticks(x)
    axC.set_xticklabels(names)
    for xi, m, col in zip(x, [x2[2] for x2 in trip], tier_colors):
        axC.text(xi, m["bands"] * 1.18,
                 f"{m['bands']} 级\n$\\Delta E$ {m['max_step_de']:.2f}",
                 ha="center", va="bottom", fontsize=9, color=INK2)
    axC.axhline(baseline_bands, color=MUTED, lw=1.2, ls="--")
    axC.annotate(f"8 位全色域参照 {baseline_bands} 级",
                 xy=(2.42, baseline_bands), xytext=(0.05, baseline_bands * 1.7),
                 color=INK2, fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9))
    axC.set_ylabel(f"b={chosen} 往返后剩余台阶数（$\\log_2$ 轴）")
    axC.set_title(f"落点检验：最窄档仅 {band_counts[0]} 级，台阶仍明显",
                  loc="left", fontsize=12, color=INK, pad=8)
    axC.grid(True, axis="y", color=GRID, lw=0.8)
    axC.set_ylim(8, max(band_counts) * 2.4)

    fig.suptitle("顶点色量化：台阶—位数—体积的落点分析（CIEDE2000 / sRGB D65）",
                 fontsize=14, color=INK, y=0.985)
    path = os.path.join(OUT, "quant_landing.png")
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    print(f"\nfigure -> {path}")


if __name__ == "__main__":
    main()
