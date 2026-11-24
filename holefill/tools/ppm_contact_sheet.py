#!/usr/bin/env python3
"""Convert holefill PPM renders to PNG and build labeled contact sheets.

Usage: python3 ppm_contact_sheet.py <dir-with-ppm>
Writes <dir>/*.png plus sheet_correct.png / sheet_buggy.png.
"""
import sys
from pathlib import Path

from PIL import Image, ImageDraw


def load_ppm(path: Path) -> Image.Image:
    return Image.open(path).convert("RGB")


def labeled(img: Image.Image, title: str, subtitle: str) -> Image.Image:
    pad = 8
    bar = 44
    w = img.width + pad * 2
    h = img.height + bar + pad * 2
    canvas = Image.new("RGB", (w, h), (24, 26, 31))
    canvas.paste(img, (pad, bar + pad))
    d = ImageDraw.Draw(canvas)
    d.text((pad, 6), title, fill=(240, 240, 240))
    d.text((pad, 24), subtitle, fill=(170, 180, 195))
    return canvas


def sheet(directory: Path, prefix: str, out_name: str, subs: dict) -> None:
    panels = []
    for kind in ("up", "side", "down"):
        p = directory / f"{prefix}_{kind}.ppm"
        if not p.exists():
            continue
        panels.append(labeled(load_ppm(p), f"{prefix} / hole {kind}", subs[kind]))
    if not panels:
        return
    gap = 10
    w = sum(p.width for p in panels) + gap * (len(panels) + 1)
    h = max(p.height for p in panels) + gap * 2
    out = Image.new("RGB", (w, h), (12, 13, 16))
    x = gap
    for p in panels:
        out.paste(p, (x, gap))
        x += p.width + gap
    out.save(directory / out_name)
    print("wrote", directory / out_name)


def main() -> None:
    directory = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    for ppm in directory.glob("*.ppm"):
        png = ppm.with_suffix(".png")
        load_ppm(ppm).save(png)
        print("wrote", png)

    good_subs = {
        "up": "dir 1.00  out 7/7  black 0  A/Aproj 1.000",
        "side": "dir 1.00  out 7/7  black 0  A/Aproj 1.000",
        "down": "dir 1.00  out 7/7  black 0  A/Aproj 1.000",
    }
    bug_subs = {
        "up": "fixed +Z agrees: black 0",
        "side": "n.ref=0 no signal: black 1804/1804",
        "down": "n.ref=-1 reversed: black 1804/1804, ratio 1.000",
    }
    sheet(directory, "correct", "sheet_correct.png", good_subs)
    sheet(directory, "buggy", "sheet_buggy.png", bug_subs)


if __name__ == "__main__":
    main()
