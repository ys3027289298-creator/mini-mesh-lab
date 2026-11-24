#!/usr/bin/env python3
"""
Size-aware normal smoothing threshold experiment for mini-mesh-lab.

The model is the same folded curved rock patch at three triangulation densities.
A fixed scanner/export position uncertainty SIGMA_Z is present at every density,
so normal error grows as faces become smaller.

Run from the repository root:
    python3 normal/threshold_lab.py
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "normal-threshold"
OUT.mkdir(parents=True, exist_ok=True)

TAU_GLOBAL_DEG = 20.0
SIGMA_Z = 0.006
CURVATURE = 0.5          # k in z = 1/2 k x^2
CREASE_SLOPE = 0.344     # |x| term; nominal normal jump about 2 atan(0.344)

# Constant linear refinement factor of three between tiers.
TIERS = [
    {"key": "large", "label": "large faces", "nx": 6, "ny": 4, "seed": 300},
    {"key": "medium", "label": "medium faces", "nx": 18, "ny": 12, "seed": 301},
    {"key": "small", "label": "small faces", "nx": 54, "ny": 36, "seed": 302},
]

# Palette: red = hard classification, amber = preserved crease that was smoothed.
C_HARD = "#c93a32"
C_FALSE = "#e69f00"
C_MESH = "#c8d0dd"
C_INK = "#20242a"
C_MUTED = "#6b7280"
C_GOOD = "#1b9e77"
C_RANGE = "#4c78a8"


def build_mesh(nx: int, ny: int, seed: int | None):
    """Build a curved patch with one exact x=0 crease and fixed z noise."""
    rng = np.random.default_rng(seed) if seed is not None else None
    xs = np.linspace(-1.0, 1.0, nx + 1)
    ys = np.linspace(-0.75, 0.75, ny + 1)
    x, y = np.meshgrid(xs, ys, indexing="ij")

    # Smooth quadratic curvature + a true |x| fold + small deterministic y form.
    z = 0.5 * CURVATURE * x * x + CREASE_SLOPE * np.abs(x)
    z += 0.03 * np.sin(3.0 * y) * x * x + 0.05 * y * y
    if rng is not None:
        z = z + rng.normal(0.0, SIGMA_Z, z.shape)

    verts = np.stack([x, y, z], axis=-1).reshape(-1, 3)

    def vid(i: int, j: int) -> int:
        return i * (ny + 1) + j

    faces: list[tuple[int, int, int]] = []
    for i in range(nx):
        for j in range(ny):
            a = vid(i, j)
            b = vid(i + 1, j)
            c = vid(i + 1, j + 1)
            d = vid(i, j + 1)
            # Same diagonal orientation in every quad.
            faces.append((a, b, c))
            faces.append((a, c, d))
    faces_arr = np.asarray(faces, dtype=np.int64)

    # Ground-truth crease: shared edges on x=0 running in y.
    i0 = nx // 2
    crease_edges = {
        tuple(sorted((vid(i0, j), vid(i0, j + 1)))) for j in range(ny)
    }

    normals = []
    for f in faces_arr:
        a, b, c = verts[f]
        n = np.cross(b - a, c - a)
        norm = np.linalg.norm(n)
        if norm == 0.0:
            raise RuntimeError("degenerate face in threshold experiment")
        normals.append(n / norm)
    normals_arr = np.asarray(normals)

    edge_faces: dict[tuple[int, int], list[int]] = {}
    for fi, f in enumerate(faces_arr):
        for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            e = tuple(sorted((int(a), int(b))))
            edge_faces.setdefault(e, []).append(fi)

    edges = []
    angles = []
    is_crease = []
    for e, fs in edge_faces.items():
        if len(fs) != 2:
            continue  # boundary edge: no adjacent-normal smoothing decision
        cosine = float(np.clip(np.dot(normals_arr[fs[0]], normals_arr[fs[1]]), -1, 1))
        theta = math.degrees(math.acos(cosine))
        edges.append(e)
        angles.append(theta)
        is_crease.append(e in crease_edges)

    edges_arr = np.asarray(edges, dtype=np.int64)
    angles_arr = np.asarray(angles)
    crease_mask = np.asarray(is_crease, dtype=bool)

    edge_lengths = np.linalg.norm(verts[edges_arr[:, 0]] - verts[edges_arr[:, 1]], axis=1)
    h = float(np.mean(edge_lengths))

    # Trace of covariance of all unit face normals.
    mean_n = normals_arr.mean(axis=0)
    global_normal_variance = float(
        np.mean(np.sum((normals_arr - mean_n) ** 2, axis=1))
    )

    # Local pair variance V = ||n1-n2||^2 / 4 = sin^2(theta/2), on smooth edges.
    smooth_theta = np.radians(angles_arr[~crease_mask])
    local_normal_variance = float(np.mean(np.sin(smooth_theta / 2.0) ** 2))

    hard_mask = angles_arr >= TAU_GLOBAL_DEG
    false_smoothed_crease_mask = crease_mask & ~hard_mask

    smooth_max = float(angles_arr[~crease_mask].max())
    crease_min = float(angles_arr[crease_mask].min())
    grouping_tau = smooth_max + 1e-9
    grouping_horn_false_creases = int(
        np.sum(crease_mask & (angles_arr < grouping_tau))
    )
    protection_horn_mask = angles_arr >= crease_min
    protection_horn_hard_edges = int(protection_horn_mask.sum())
    protection_horn_hard_smooth_edges = int(
        np.sum(protection_horn_mask & ~crease_mask)
    )

    return {
        "tier": None,
        "nx": nx,
        "ny": ny,
        "h_mean_edge": h,
        "vertices": int(len(verts)),
        "faces": int(len(faces_arr)),
        "interior_edges": int(len(edges_arr)),
        "crease_edges_truth": int(crease_mask.sum()),
        "hard_edges": int(hard_mask.sum()),
        "false_smoothed_creases": int(false_smoothed_crease_mask.sum()),
        "normal_variance": local_normal_variance,
        "global_face_normal_variance": global_normal_variance,
        "smooth_angle_min_deg": float(angles_arr[~crease_mask].min()),
        "smooth_angle_max_deg": float(angles_arr[~crease_mask].max()),
        "smooth_angle_mean_deg": float(angles_arr[~crease_mask].mean()),
        "smooth_angle_p99_deg": float(np.percentile(angles_arr[~crease_mask], 99)),
        "crease_angle_min_deg": crease_min,
        "crease_angle_max_deg": float(angles_arr[crease_mask].max()),
        "feasible_zero_error": bool(smooth_max < crease_min),
        "grouping_horn_false_creases": grouping_horn_false_creases,
        "protection_horn_hard_edges": protection_horn_hard_edges,
        "protection_horn_hard_smooth_edges": protection_horn_hard_smooth_edges,
        "_verts": verts,
        "_faces": faces_arr,
        "_normals": normals_arr,
        "_edges": edges_arr,
        "_angles": angles_arr,
        "_crease_mask": crease_mask,
        "_hard_mask": hard_mask,
        "_false_mask": false_smoothed_crease_mask,
    }


def render_tier(data: dict, path: Path) -> None:
    verts = data["_verts"]
    faces = data["_faces"]
    normals = data["_normals"]
    edges = data["_edges"]
    hard_mask = data["_hard_mask"]
    false_mask = data["_false_mask"]

    fig = plt.figure(figsize=(12.8, 5.3), dpi=170)
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax2d = fig.add_subplot(1, 2, 2)

    fig.suptitle(
        f"{data['tier']} faces: h={data['h_mean_edge']:.3f} | "
        f"hard edges={data['hard_edges']} | "
        f"false-smoothed creases={data['false_smoothed_creases']} | "
        f"V_N={data['normal_variance']:.6f}",
        x=0.035,
        y=0.97,
        ha="left",
        color=C_INK,
        fontsize=12,
        fontweight="bold",
    )

    ax3d.plot_trisurf(
        verts[:, 0],
        verts[:, 1],
        verts[:, 2],
        triangles=faces,
        color=C_MESH,
        edgecolor="none",
        linewidth=0.0,
        antialiased=True,
        shade=True,
    )

    lw_hard = 1.8 if data["nx"] <= 6 else (1.05 if data["nx"] <= 18 else 0.55)
    hard_lc3 = Line3DCollection(
        verts[edges[hard_mask]],
        colors=C_HARD,
        linewidths=lw_hard,
        alpha=0.95,
    )
    hard_lc3.set_zorder(5)
    ax3d.add_collection3d(hard_lc3)

    false_lw = 2.5 if data["nx"] <= 18 else 1.8
    if false_mask.any():
        false_lc3 = Line3DCollection(
            verts[edges[false_mask]],
            colors=C_FALSE,
            linewidths=false_lw,
            alpha=1.0,
        )
        false_lc3.set_zorder(10)
        ax3d.add_collection3d(false_lc3)

    ax3d.set_xlim(-1.05, 1.05)
    ax3d.set_ylim(-0.80, 0.80)
    ax3d.set_zlim(-0.08, 0.82)
    ax3d.set_box_aspect((2.1, 1.5, 0.72))
    ax3d.view_init(elev=25, azim=-58)
    ax3d.set_axis_off()
    ax3d.set_title("shaded oblique view", fontsize=10, color=C_MUTED, pad=0)

    # Top view with explicit face shading, so edge density remains readable.
    light = np.array([-0.35, -0.45, 0.82])
    light /= np.linalg.norm(light)
    diffuse = np.clip(normals @ light, 0.0, 1.0)
    base = np.array([0.76, 0.80, 0.87])
    face_colors = base[None, :] * (0.60 + 0.40 * diffuse[:, None])
    polygons = verts[faces][:, :, :2]
    face_collection = PolyCollection(
        polygons,
        facecolors=face_colors,
        edgecolors="none",
        zorder=1,
    )
    ax2d.add_collection(face_collection)

    hard_lc2 = LineCollection(
        verts[edges[hard_mask], :2],
        colors=C_HARD,
        linewidths=lw_hard,
        alpha=0.9,
    )
    hard_lc2.set_zorder(4)
    ax2d.add_collection(hard_lc2)
    if false_mask.any():
        false_lc2 = LineCollection(
            verts[edges[false_mask], :2],
            colors=C_FALSE,
            linewidths=false_lw,
        )
        false_lc2.set_zorder(6)
        ax2d.add_collection(false_lc2)

    ax2d.set_xlim(-1.03, 1.03)
    ax2d.set_ylim(-0.78, 0.78)
    ax2d.set_aspect("equal")
    ax2d.set_xticks([])
    ax2d.set_yticks([])
    ax2d.set_title("top view: hard-edge classification", fontsize=10, color=C_MUTED)
    for spine in ax2d.spines.values():
        spine.set_color("#d1d5db")

    legend_handles = [
        Line2D([0], [0], color=C_HARD, lw=2.2, label=f"hard at τ={TAU_GLOBAL_DEG:g}°"),
        Line2D([0], [0], color=C_FALSE, lw=2.8, label="true crease smoothed"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.925),
        frameon=True,
        framealpha=0.94,
        fontsize=9,
    )
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.03, top=0.80, wspace=0.02)
    fig.savefig(path, pad_inches=0.12, facecolor="white")
    plt.close(fig)


def render_interval_chart(rows: list[dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 3.9), dpi=170)
    y = np.arange(len(rows))[::-1]

    for yi, r in zip(y, rows):
        lower = r["smooth_angle_max_deg"]
        upper = r["crease_angle_min_deg"]
        if lower < upper:
            ax.plot([lower, upper], [yi, yi], color=C_GOOD, lw=8, solid_capstyle="butt")
            ax.text(
                (lower + upper) / 2,
                yi + 0.23,
                f"usable ({lower:.1f}°, {upper:.1f}°)",
                ha="center",
                va="bottom",
                fontsize=9,
                color=C_INK,
            )
        else:
            # Show the contradictory required rays rather than a fake interval.
            ax.plot([lower, lower + 18], [yi, yi], color=C_HARD, lw=8, solid_capstyle="butt")
            ax.plot([0, upper], [yi, yi], color=C_FALSE, lw=8, solid_capstyle="butt")
            ax.text(
                lower + 1,
                yi + 0.24,
                f"need τ>{lower:.1f}°",
                ha="left",
                va="bottom",
                fontsize=9,
                color=C_HARD,
            )
            ax.text(
                upper + 2.0,
                yi - 0.28,
                f"need τ<{upper:.1f}°",
                ha="left",
                va="top",
                fontsize=9,
                color="#9a6500",
            )

    ax.axvline(
        TAU_GLOBAL_DEG,
        color=C_INK,
        ls="--",
        lw=1.4,
        label=f"global τ={TAU_GLOBAL_DEG:g}°",
    )
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r['tier']}\nh={r['h_mean_edge']:.3f}" for r in rows], fontsize=9
    )
    ax.set_xlabel("adjacent face-normal angle θ (degrees)")
    ax.set_xlim(0, 82)
    ax.set_ylim(-0.7, len(rows) - 0.1)
    ax.grid(axis="x", color="#e5e7eb", lw=0.8)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    rows = []
    for spec in TIERS:
        data = build_mesh(spec["nx"], spec["ny"], spec["seed"])
        data["tier"] = spec["key"]
        render_tier(data, OUT / f"render_{spec['key']}.png")
        rows.append(data)

    render_interval_chart(rows, OUT / "threshold_intervals.png")

    table_rows = []
    for r in rows:
        table_rows.append(
            {
                "tier": r["tier"],
                "h_mean_edge": r["h_mean_edge"],
                "hard_edges": r["hard_edges"],
                "false_smoothed_creases": r["false_smoothed_creases"],
                "normal_variance": r["normal_variance"],
            }
        )

    with (OUT / "nine_numbers.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(table_rows[0].keys()))
        writer.writeheader()
        writer.writerows(table_rows)

    serializable = []
    for r in rows:
        serializable.append({k: v for k, v in r.items() if not k.startswith("_")})
    payload = {
        "definition": {
            "threshold_degrees": TAU_GLOBAL_DEG,
            "hard_when": "theta >= tau",
            "normal_variance": "mean over smooth shared edges of ||n1-n2||^2/4 = sin^2(theta/2)",
            "position_noise_sigma_z": SIGMA_Z,
            "zero_error_interval": "(max smooth-edge angle, min crease-edge angle)",
        },
        "tiers": serializable,
    }
    (OUT / "metrics.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("tier      h      hard false  V_normal  L_smooth  U_crease feasible")
    for r in rows:
        print(
            f"{r['tier']:<6} {r['h_mean_edge']:.4f} "
            f"{r['hard_edges']:5d} {r['false_smoothed_creases']:5d} "
            f"{r['normal_variance']:.6f}  {r['smooth_angle_max_deg']:8.3f} "
            f"{r['crease_angle_min_deg']:8.3f} {r['feasible_zero_error']}"
        )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
