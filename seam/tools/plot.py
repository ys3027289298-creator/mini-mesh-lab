#!/usr/bin/env python3
"""Pure-stdlib SVG plot of the luminance sampling lines (no matplotlib)."""
import csv
import sys


def read(path):
    cases = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            cases.setdefault(row["case"], []).append(
                (float(row["s"]), float(row["L_flat"]),
                 float(row["L_area"]), float(row["L_fixed"])))
    return cases


def polyline(points, sx, sy, col):
    pts = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in points)
    return f'<polyline fill="none" stroke="{col}" stroke-width="1.6" points="{pts}"/>'


def panel(case, rows, x0, y0, w, h, title):
    s = [r[0] for r in rows]
    smin, smax = min(s), max(s)
    seam1, seam2 = 2.0, smax - 2.0          # flat length is 2 on each side
    sx = lambda x: x0 + (x - smin) / (smax - smin) * w
    sy = lambda y: y0 + h - y * h
    out = [f'<text x="{x0}" y="{y0 - 8}" font-size="14" font-family="sans-serif">{title}</text>']
    for sxm in (seam1, seam2):
        xx = sx(sxm)
        out.append(f'<line x1="{xx:.1f}" y1="{y0}" x2="{xx:.1f}" y2="{y0 + h}" '
                   f'stroke="#888" stroke-dasharray="4,3" stroke-width="1"/>')
    out.append(f'<text x="{sx(seam1):.1f}" y="{y0 + h - 6}" font-size="10" '
               f'font-family="sans-serif" fill="#666">seam</text>')
    series = [("L_flat", 1, "#c0392b"), ("L_area (current)", 2, "#e67e22"),
              ("L_fixed", 3, "#1e8449")]
    for name, idx, col in series:
        out.append(polyline([(r[0], r[idx]) for r in rows], sx, sy, col))
    for gy in range(0, 6):
        yy = y0 + h - gy / 5 * h
        out.append(f'<line x1="{x0}" y1="{yy:.1f}" x2="{x0 + w}" y2="{yy:.1f}" '
                   f'stroke="#eee" stroke-width="0.8"/>')
        out.append(f'<text x="{x0 - 30}" y="{yy + 4:.1f}" font-size="10" '
                   f'font-family="sans-serif" fill="#666">{gy / 5:.1f}</text>')
    out.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#333"/>')
    return "\n".join(out)


def main():
    cases = read(sys.argv[1])
    w, h = 720, 210
    panels = []
    legend_y = 26
    legend = ['<text x="40" y="20" font-size="12" font-family="sans-serif">luminance along sampling line '
              '(dashed grey = bevel seams)</text>']
    for i, (name, col) in enumerate([("L_flat", "#c0392b"),
                                     ("L_area (current)", "#e67e22"),
                                     ("L_fixed", "#1e8449")]):
        legend.append(f'<line x1="{300 + i * 150}" y1="16" x2="{324 + i * 150}" y2="16" '
                      f'stroke="{col}" stroke-width="2"/>')
        legend.append(f'<text x="{330 + i * 150}" y="20" font-size="12" '
                      f'font-family="sans-serif">{name}</text>')
    for j, (tag, title) in enumerate([("uniform8", "uniform 8-segment bevel"),
                                      ("nonuniform8", "non-uniform 8-segment bevel")]):
        panels.append(panel(tag, cases[tag], 70, 90 + j * (h + 60), w, h, title))
    total_h = 90 + 2 * (h + 60)
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="860" height="{total_h}" '
           f'viewBox="0 0 860 {total_h}">\n' + "\n".join(legend + panels) + "\n</svg>\n")
    with open(sys.argv[2], "w") as f:
        f.write(svg)
    print(f"wrote {sys.argv[2]}")


if __name__ == "__main__":
    main()
