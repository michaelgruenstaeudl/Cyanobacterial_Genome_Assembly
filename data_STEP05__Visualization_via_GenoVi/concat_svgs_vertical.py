#!/usr/bin/env python3
"""Concatenate multiple SVG files vertically into a single SVG."""

from __future__ import annotations

import argparse
import copy
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from cairosvg import svg2png

SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Concatenate SVGs vertically (top to bottom)."
    )
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        required=True,
        type=Path,
        help="Input SVG files in top-to-bottom order.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output basename or SVG path (e.g., result or result.svg).",
    )
    return parser.parse_args()


def _parse_length(value: str | None) -> float | None:
    if not value:
        return None
    m = re.match(r"^\s*([+-]?(?:\d+\.?\d*|\d*\.\d+))", value)
    if not m:
        return None
    return float(m.group(1))


def _get_dimensions(root: ET.Element, src: Path) -> tuple[float, float]:
    viewbox = root.attrib.get("viewBox")
    if viewbox:
        parts = viewbox.replace(",", " ").split()
        if len(parts) == 4:
            return float(parts[2]), float(parts[3])

    width = _parse_length(root.attrib.get("width"))
    height = _parse_length(root.attrib.get("height"))
    if width is None or height is None:
        raise ValueError(f"Could not determine SVG dimensions for {src}")
    return width, height


def _load_svg(path: Path) -> tuple[ET.Element, float, float]:
    tree = ET.parse(path)
    root = tree.getroot()
    width, height = _get_dimensions(root, path)
    return root, width, height


def main() -> None:
    args = parse_args()
    base_out = args.output
    if base_out.suffix.lower() in (".svg", ".png"):
        base_out = base_out.with_suffix("")
    svg_out = base_out.with_suffix(".svg")
    png_out = svg_out.with_suffix(".png")

    svgs = []
    for input_path in args.inputs:
        if not input_path.exists():
            raise FileNotFoundError(f"Input SVG not found: {input_path}")
        svgs.append(_load_svg(input_path))

    out_width = max(w for _, w, _ in svgs)
    scales = [out_width / w if w else 1.0 for _, w, _ in svgs]
    scaled_heights = [h * s for (_, _, h), s in zip(svgs, scales)]
    out_height = sum(scaled_heights)

    out_root = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "version": "1.1",
            "width": str(out_width),
            "height": str(out_height),
            "viewBox": f"0 0 {out_width} {out_height}",
        },
    )

    y_offset = 0.0
    for (root, width, height), scale in zip(svgs, scales):
        group = ET.SubElement(
            out_root,
            f"{{{SVG_NS}}}g",
            {"transform": f"translate(0,{y_offset:.6f}) scale({scale:.6f})"},
        )
        for child in list(root):
            group.append(copy.deepcopy(child))
        y_offset += height * scale

    ET.ElementTree(out_root).write(svg_out, encoding="utf-8", xml_declaration=True)
    svg2png(url=str(svg_out), write_to=str(png_out))
    print(f"Saved concatenated SVG to: {svg_out}")
    print(f"Saved concatenated PNG to: {png_out}")


if __name__ == "__main__":
    main()
