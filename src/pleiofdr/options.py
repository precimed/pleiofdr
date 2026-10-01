"""Analysis options (port of pleioOpt.m)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

STATTYPES = ("condfdr", "conjfdr")
REPEATS = ("default", "maxout", "none")


def matlab_linspace(d1: float, d2: float, n: int) -> np.ndarray:
    """MATLAB linspace: d1 + (0:n-1)*(d2-d1)/(n-1), multiplying before dividing.

    numpy multiplies by a precomputed step instead, which moves some grid points by one ulp
    (e.g. 29.99 vs 29.990000000000002) and with them the histogram bin edges.
    """
    n = int(n)
    if n < 2:
        return np.array([float(d2)] * max(n, 0))
    y = d1 + np.arange(n) * (d2 - d1) / (n - 1)
    y[0], y[-1] = d1, d2
    return y


def _default_colorlist() -> np.ndarray:
    return 0.8 * np.array(
        [[1.00, 0, 0], [1.00, 0.5, 0], [0, 0.75, 0.75], [0, 0.50, 0], [0.75, 0, 0.75], [0, 0, 1.00]]
    )


@dataclass
class Options:
    stattype: str = "condfdr"
    fdrthresh: float = 0.05
    pthresh: float = 5e-8
    # 0-based inclusive (first, last) SNP index pairs; pleioOpt stores 1-based pairs
    exclude_from_fit: list[tuple[int, int]] = field(default_factory=list)
    mafthresh: float = 0.005
    t1_low: float = 0
    t1_up: float = 30
    t1_nbreaks: int = 3001
    t2_low: float = 0
    t2_up: float = 3
    t2_nbreaks: int = 31
    qq_low: float = 0
    qq_up: float = 3
    qq_nbreaks: int = 4
    onscreen: bool = False
    randprune: bool = True
    randprune_n: int = 100
    randprune_file: str = ""
    randprune_repeats: str = "default"
    reset_pruneidx: bool = True
    smf: float = 1e2
    thin: int = 10
    perform_gc: bool = True
    randprune_gc: bool = False
    smooth_lookup: bool = True
    adjust_lookup: bool = False
    show_ci: bool = False
    use_standard_gc: bool = False
    exclude_ambiguous_snps: bool = False
    dummy_zscore: bool = False
    exclude_from_fit_and_discovery: bool = False
    outputdir: str = ""
    correct_for_sample_overlap: bool = False
    fishercomb: bool = True
    manh_fontsize_genenames: float = 20
    manh_fontsize_legends: float = 18
    manh_fontsize_axes: float = 18
    manh_legend: str = "NorthEast"
    manh_plot: bool = False
    manh_ymargin: float = 1
    manh_yspace: float = 0.75
    manh_image_size: tuple[float, float] | None = (1920, 1440)
    manh_colorlist: np.ndarray = field(default_factory=_default_colorlist)
    seed: int | None = None

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if self.stattype.lower() not in STATTYPES:
            raise ValueError(f"Not-implemented FDR type {self.stattype}")
        self.stattype = self.stattype.lower()
        if self.randprune_repeats not in REPEATS:
            raise ValueError(f"randprune_repeats must be one of {REPEATS}")
        for first, last in self.exclude_from_fit:
            if first >= last:
                raise ValueError("exclude region is empty")
            if first < 0 or last < 0:
                raise ValueError("Invalid exclude region: negative values")

    @property
    def qqbreaks(self) -> np.ndarray:
        return matlab_linspace(self.qq_low, self.qq_up, self.qq_nbreaks)

    @property
    def t1breaks(self) -> np.ndarray:
        return matlab_linspace(self.t1_low, self.t1_up, self.t1_nbreaks)

    @property
    def t2breaks(self) -> np.ndarray:
        return matlab_linspace(self.t2_low, self.t2_up, self.t2_nbreaks)

    @property
    def hv(self) -> np.ndarray:
        """Thinned trait-1 grid, opts.t1breaks(1:opts.thin:end)."""
        return self.t1breaks[:: self.thin]
