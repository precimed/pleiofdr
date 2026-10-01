"""Figures (ports of save_figure.m, filter_points_for_plotting.m and the plot_* functions)."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np

# MATLAB's default axes color order (R2014b and later)
MATLAB_COLORS = np.array(
    [
        [0.0000, 0.4470, 0.7410],
        [0.8500, 0.3250, 0.0980],
        [0.9290, 0.6940, 0.1250],
        [0.4940, 0.1840, 0.5560],
        [0.4660, 0.6740, 0.1880],
        [0.3010, 0.7450, 0.9330],
        [0.6350, 0.0780, 0.1840],
    ]
)


def use_backend(onscreen: bool) -> None:
    """Non-interactive rendering unless figures should be shown on screen."""
    if not onscreen:
        matplotlib.use("Agg", force=True)


_TEX_ESCAPES = {
    "\\": r"\backslash ",
    "_": r"\_",
    "%": r"\%",
    "$": r"\$",
    "{": r"\{",
    "}": r"\}",
    "^": r"\^{}",
    "~": r"\sim ",
    " ": r"\ ",
    "#": r"\text{\#}",
}


def tex(name: str) -> str:
    """Escape a trait name for use inside matplotlib mathtext."""
    return "".join(_TEX_ESCAPES.get(ch, ch) for ch in name)


def save_figure(fig, fmt: str, filename: str | Path) -> Path | None:
    """Write fig as <filename>.<fmt>; MATLAB's .fig format has no equivalent and is skipped."""
    if fig is None or fmt == "fig":
        return None
    stem = str(filename)
    if stem.endswith(f".{fmt}"):
        stem = stem[: -len(fmt) - 1]
    outfn = Path(f"{stem}.{fmt}")
    print(f"writing {outfn}")
    fig.savefig(outfn, format=fmt)
    return outfn


def filter_points(
    x: np.ndarray, y: np.ndarray, image_size: tuple[float, float] | None
) -> tuple[np.ndarray, np.ndarray]:
    """Keep the first point in each pixel of an image_size raster; drops NaN points."""
    if image_size is None or any(np.isnan(image_size)):
        return x, y
    if len(x) != len(y):
        raise ValueError("x and y length mismatch")
    keep = ~(np.isnan(x) | np.isnan(y))
    x, y = x[keep], y[keep]
    if x.size == 0:
        return x, y
    width, height = int(image_size[0]), int(image_size[1])

    def pixel(v: np.ndarray, n: int) -> np.ndarray:
        span = (v.max() - v.min()) / n
        if span == 0:
            return np.ones(v.size, dtype=np.int64)  # int32(NaN) is 0 in MATLAB
        return 1 + np.floor((v - v.min()) / span + 0.5).astype(np.int64)  # int32() rounds

    key = pixel(x, width) * (height + 2) + pixel(y, height)
    _, first = np.unique(key, return_index=True)
    first.sort()
    return x[first], y[first]


def subplot_grid(n: int) -> tuple[int, int]:
    cols = int(np.ceil(np.sqrt(n)))
    return int(np.ceil(n / cols)), cols
