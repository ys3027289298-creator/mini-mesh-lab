#!/usr/bin/env python3
"""Build the seam-tolerance figures from the experiment CSVs.

Reads sweep.csv / summary.csv / scale.csv from the artifacts directory given
in argv[1] and writes PNG figures next to them.  Palette is the validated
dataviz reference set (categorical slots 1-4, light surface)."""
import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
RED = "#e34948"
VIOLET = "#6250d6"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "system-ui", "Segoe UI"],
    "font.size": 10,
    "text.color": INK,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK2,
    "axes.titlecolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.axisbelow": True,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

DENSITIES = ["sparse", "medium", "dense"]
DENSITY_TITLE = {"sparse": "sparse (2 seams)", "medium": "medium (7 seams)",
                 "dense": "dense (14 seams)"}


def load_sweep(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            r["k"] = float(r["k"])
            for key in ("unmerged_seam_pairs", "mismerged_crease_pairs",
                        "bright_pixels"):
                r[key] = int(r[key])
            r["mean_normal_error_deg"] = float(r["mean_normal_error_deg"])
            rows.append(r)
    return rows


def load_summary(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def tolerance_figure(rows, out):
    """The two curves that set the tolerance together:
    unmerged seams fall as k grows, mismerged creases rise."""
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=False)
    for ax, dens in zip(axes, DENSITIES):
        sub = [r for r in rows if r["density"] == dens and r["mode"] == "position"]
        sub.sort(key=lambda r: r["k"])
        ks = [r["k"] for r in sub]
        unmerged = [r["unmerged_seam_pairs"] for r in sub]
        mismerged = [r["mismerged_crease_pairs"] for r in sub]
        ymax = max(max(unmerged), max(mismerged)) * 1.18
        ax.semilogx(ks, unmerged, color=BLUE, lw=2, label="unmerged seam pairs")
        ax.semilogx(ks, mismerged, color=ORANGE, lw=2, label="mismerged crease pairs")
        k_seam = max(r["k"] for r in sub if r["unmerged_seam_pairs"] > 0)
        k_crease = min(r["k"] for r in sub if r["mismerged_crease_pairs"] > 0)
        k_star = (k_seam * k_crease) ** 0.5
        marks = ((k_seam, BLUE, "k_seam", 0.98),
                 (k_star, VIOLET, "k*", 0.86),
                 (k_crease, ORANGE, "k_crease", 0.74))
        for x, c, txt, frac in marks:
            ax.axvline(x, color=c, lw=1.2, ls="--", alpha=0.8)
            ax.annotate(txt, xy=(x, ymax * frac), fontsize=8, color=c,
                        ha="center", va="top", rotation=90,
                        xytext=(0, -2), textcoords="offset points")
        ax.set_title(DENSITY_TITLE[dens], fontsize=11)
        ax.set_xlabel("tolerance k = eps / L")
        ax.set_ylim(-4, ymax)
        # direct end labels as secondary encoding for the CVD warn band
        ax.annotate("unmerged", xy=(ks[-1], unmerged[-1]), color=BLUE,
                    fontsize=8, va="center", xytext=(5, 8),
                    textcoords="offset points")
        ax.annotate("mismerged", xy=(ks[-1], mismerged[-1]), color=ORANGE,
                    fontsize=8, va="center", xytext=(5, -10),
                    textcoords="offset points")
    axes[0].set_ylabel("vertex-pair count")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Tolerance vs patch size: the two competing failure curves",
                 fontsize=13, y=1.02)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def error_figure(rows, out):
    """Normal error and bright pixels vs k, position vs guarded weld."""
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2), sharex=True)
    for col, dens in enumerate(DENSITIES):
        for mode, color in (("position", BLUE), ("guarded", AQUA)):
            sub = sorted((r for r in rows if r["density"] == dens
                          and r["mode"] == mode), key=lambda r: r["k"])
            ks = [r["k"] for r in sub]
            axes[0][col].semilogx(ks, [r["mean_normal_error_deg"] for r in sub],
                                  color=color, lw=2, label=mode)
            axes[1][col].semilogx(ks, [r["bright_pixels"] for r in sub],
                                  color=color, lw=2, label=mode)
        axes[0][col].set_title(DENSITY_TITLE[dens], fontsize=11)
        axes[1][col].set_xlabel("tolerance k = eps / L")
    axes[0][0].set_ylabel("mean normal error (deg)")
    axes[1][0].set_ylabel("bright seam pixels")
    axes[0][0].legend(frameon=False, fontsize=9)
    fig.suptitle("Normal error and visible seam lines across the tolerance sweep",
                 fontsize=13, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def metrics_figure(summary, out):
    """The headline 3 densities x 4 metrics, two normal-field variants."""
    metrics = [
        ("unmerged_vertices", "unmerged seam vertices"),
        ("mismerged_crease_pairs", "mismerged crease pairs"),
        ("bright_pixels", "bright seam pixels"),
        ("mean_normal_error_deg", "mean normal error (deg)"),
    ]
    variants = [("position_match_only", BLUE, "position match only"),
                ("channel_open_position", AQUA, "channel opened (shared normal)")]
    fig, axes = plt.subplots(1, 4, figsize=(14, 4.2))
    width = 0.36
    x = range(len(DENSITIES))
    for ax, (key, title) in zip(axes, metrics):
        for vi, (variant, color, label) in enumerate(variants):
            vals = []
            for dens in DENSITIES:
                r = next(r for r in summary
                         if r["density"] == dens and r["scale"] == "1.00"
                         and r["variant"] == variant)
                vals.append(float(r[key]))
            offset = -width / 2 if vi == 0 else width / 2
            bars = ax.bar([i + offset for i in x], vals, width=width,
                          color=color, label=label)
            for b, v in zip(bars, vals):
                ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                        f"{v:g}", ha="center", va="bottom", fontsize=8,
                        color=INK2)
        ax.set_xticks(list(x))
        ax.set_xticklabels(["2 seams", "7 seams", "14 seams"], fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.set_ylim(bottom=0)
        # keep all-zero panels readable instead of a symmetric auto range
        top = ax.get_ylim()[1]
        if top < 1.0:
            ax.set_ylim(0, 1.0)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Three seam densities, four metrics (unit scale, k = k*)",
                 fontsize=13, y=1.02)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def channel_figure(summary, out):
    """If the densest model still shows seam lines, compare UV channels
    against the number of distinct normals per weld group."""
    fig, ax = plt.subplots(figsize=(9, 4.6))
    width = 0.26
    x = range(len(DENSITIES))
    series = [
        ("uv_channels", MUTED, "UV islands (topology)"),
        ("normal_values", BLUE, "distinct normals / weld group: position match only"),
        ("normal_values_open", AQUA, "distinct normals / weld group: channel opened"),
    ]
    for si, (key, color, label) in enumerate(series):
        vals = []
        for dens in DENSITIES:
            if key == "normal_values_open":
                variant = "channel_open_position"
                col = "normal_values"
            elif key == "normal_values":
                variant = "position_match_only"
                col = "normal_values"
            else:
                variant = "position_match_only"
                col = "uv_channels"
            r = next(r for r in summary
                     if r["density"] == dens and r["scale"] == "1.00"
                     and r["variant"] == variant)
            vals.append(int(r[col]))
        offset = (si - 1) * width
        # log scale so the handful of UV islands stays visible next to the
        # ~1000 normal channels; bars start at 1 and reach the true value
        bars = ax.bar([i + offset for i in x], [v - 1.0 for v in vals],
                      width=width, bottom=1.0, color=color, label=label)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() * 1.08,
                    str(v), ha="center", va="bottom", fontsize=8, color=INK2)
    ax.set_yscale("log")
    ax.set_ylim(1, 4000)
    ax.set_xticks(list(x))
    ax.set_xticklabels([DENSITY_TITLE[d] for d in DENSITIES])
    ax.set_ylabel("channel count (log scale)")
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("UV channels vs normal channels: position matching alone "
                 "keeps the normals split", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    art = sys.argv[1] if len(sys.argv) > 1 else "../artifacts/seam"
    rows = load_sweep(f"{art}/sweep.csv")
    summary = load_summary(f"{art}/summary.csv")
    tolerance_figure(rows, f"{art}/fig_tolerance_curves.png")
    error_figure(rows, f"{art}/fig_error_curves.png")
    metrics_figure(summary, f"{art}/fig_metrics.png")
    channel_figure(summary, f"{art}/fig_channels.png")
    print("figures written to", art)


if __name__ == "__main__":
    main()
