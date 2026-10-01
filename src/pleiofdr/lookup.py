"""Conditional FDR lookup tables (ports of lookup_table.m, conj_lookup_table.m, cond_FDR_amd.m)."""

from __future__ import annotations

import numpy as np

from .options import Options
from .smoothing import sparse_smooth_2d
from .stats import binofit


def histc_bins(x: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """0-based MATLAB histc bin of each value, -1 when it is not counted.

    Bin k holds edges[k] <= x < edges[k+1]; the last bin holds x == edges[-1]; NaN and values
    outside the edges are ignored.
    """
    k = np.searchsorted(edges, x, side="right") - 1
    last = edges.size - 1
    valid = (k >= 0) & ((k < last) | (x == edges[-1])) & ~np.isnan(x)
    return np.where(valid, k, -1)


def histc(x: np.ndarray, edges: np.ndarray) -> np.ndarray:
    k = histc_bins(x, edges)
    return np.bincount(k[k >= 0], minlength=edges.size).astype(float)


def hist3(y: np.ndarray, x: np.ndarray, edges_y: np.ndarray, edges_x: np.ndarray) -> np.ndarray:
    """MATLAB hist3([y x], 'Edges', {edges_y, edges_x})."""
    ky, kx = histc_bins(y, edges_y), histc_bins(x, edges_x)
    ok = (ky >= 0) & (kx >= 0)
    flat = ky[ok] * edges_x.size + kx[ok]
    counts = np.bincount(flat, minlength=edges_y.size * edges_x.size)
    return counts.reshape(edges_y.size, edges_x.size).astype(float)


def interp2(v: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """MATLAB interp2(V, x, y) (linear) on the 1-based grid 1:ncols x 1:nrows; NaN outside.

    Like griddedInterpolant, corners with zero weight do not contribute, so a NaN neighbour
    does not spoil a query that falls exactly on a grid line.
    """
    nr, nc = v.shape
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    inside = (x >= 1) & (x <= nc) & (y >= 1) & (y <= nr)
    xs = np.where(inside, x, 1.0)
    ys = np.where(inside, y, 1.0)
    x0 = np.clip(np.floor(xs), 1, max(nc - 1, 1)).astype(np.intp)
    y0 = np.clip(np.floor(ys), 1, max(nr - 1, 1)).astype(np.intp)
    tx, ty = xs - x0, ys - y0
    x1 = np.minimum(x0 + 1, nc)
    y1 = np.minimum(y0 + 1, nr)
    out = np.zeros(x.shape)
    for yy, xx, w in (
        (y0, x0, (1 - tx) * (1 - ty)),
        (y0, x1, tx * (1 - ty)),
        (y1, x0, (1 - tx) * ty),
        (y1, x1, tx * ty),
    ):
        out += np.where(w != 0, w * v[yy - 1, xx - 1], 0.0)
    out[~inside] = np.nan
    return out


def _logit(p: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(p / (1 - p))


def ind_look(lp1: np.ndarray, lp2: np.ndarray, opts: Options) -> np.ndarray:
    """FDR lookup table (t2_nbreaks x thinned t1 grid) from one set of SNPs."""
    t1, t2 = opts.t1breaks, opts.t2breaks
    ok = ~(np.isnan(lp1) | np.isnan(lp2))
    edges_y = np.append(t2, np.inf)  # Inf captures lp2 > max(t2breaks)
    edges_x = np.append(t1, np.inf) - 0.5 * (t1[1] - t1[0])  # half-bin shift of the legacy code
    hlp = hist3(lp2[ok], lp1[ok], edges_y, edges_x)

    hcmat12 = np.cumsum(hlp[::-1], axis=0)[::-1]  # SNPs with lp2 >= t2breaks(i)
    x = np.cumsum(hcmat12[:, ::-1], axis=1)[:, ::-1]  # ... and lp1 >= bin j
    n = np.repeat(hcmat12.sum(axis=1, keepdims=True), x.shape[1], axis=1)
    x, n = x[:-1, :-1], n[:-1, :-1]

    pest, pci = binofit(1 + x.ravel(order="F"), 1 + n.ravel(order="F"))
    pest = pest.reshape(x.shape, order="F")
    pupp = pci[:, 1].reshape(x.shape, order="F")
    plow = pci[:, 0].reshape(x.shape, order="F")

    muim = _logit(pest)
    if opts.smooth_lookup:
        with np.errstate(divide="ignore", invalid="ignore"):
            wim = (_logit(pupp) - _logit(plow)) ** -2.0
        muim = muim[:, :: opts.thin].copy()
        wim = wim[:, :: opts.thin].copy()
        bad = ~np.isfinite(muim)
        wim[bad] = 0
        muim[bad] = 0
        muim_sm = sparse_smooth_2d(muim, wim, (opts.smf, opts.smf))
    else:
        muim_sm = muim[:, :: opts.thin]
        print("WARNING!!! no smoothing")

    qim = 1.0 / (1.0 + np.exp(-muim_sm))  # observed F(p1 | p2)
    pim = np.broadcast_to(10.0 ** (-opts.hv), qim.shape)  # F0 under the null
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = pim / qim
    # MATLAB min/max ignore NaN
    return np.fmin(1.0, np.fmax(0.0, ratio))


def lookup_table(
    lp1: np.ndarray,
    lp2: np.ndarray,
    opts: Options,
    pruneidx: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    """Conditional FDR table, averaged over random-pruning iterations when pruneidx is given."""
    if pruneidx is None:
        looktable, lookcount, lookstd = ind_look(lp1, lp2, opts), None, None
    else:
        if pruneidx.shape[0] != lp1.shape[0] or pruneidx.shape[0] != lp2.shape[0]:
            raise ValueError("Random prune index does not conform to logpvec")
        shape = (opts.t2_nbreaks, opts.hv.size)
        looktable, lookcount, lookvar = np.zeros(shape), np.zeros(shape), np.zeros(shape)
        niter = int(opts.randprune_n)
        for i in range(niter):
            print(f"\rFill lookup table {i + 1:3d}/{niter:3d} ", end="", flush=True)
            keep = pruneidx[:, i]
            tmp_lp1 = np.where(keep, lp1, np.nan)
            tmp_lp2 = np.where(keep, lp2, np.nan)
            tmp_look = ind_look(tmp_lp1, tmp_lp2, opts)
            defmat = np.isfinite(tmp_look)
            lookcount += defmat
            looktable[defmat] += tmp_look[defmat]
            lookvar[defmat] += tmp_look[defmat] ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            looktable = looktable / lookcount
            # biased estimator; MATLAB returns a complex number for tiny negative variances
            lookstd = np.sqrt(np.maximum(lookvar / lookcount - looktable**2, 0.0))

    if opts.adjust_lookup:
        print("adjusting lookup... ", end="")
        for col in range(looktable.shape[1]):
            for row in range(1, looktable.shape[0]):
                looktable[row, col] = min(looktable[row, col], looktable[row - 1, col])
            if col > 0:
                looktable[:, col] = np.minimum(looktable[:, col], looktable[:, col - 1])
    return looktable, lookcount, lookstd


def conj_lookup_table(
    lp1: np.ndarray, lp2: np.ndarray, opts: Options, pruneidx: np.ndarray | None = None
) -> np.ndarray:
    """conj_lookup_table.m: note both tables are built as lp1 | lp2, as in the MATLAB code."""
    look12 = lookup_table(lp1, lp2, opts, pruneidx)[0]
    look21 = lookup_table(lp1, lp2, opts, pruneidx)[0]
    return combine_conj_lookup(look12, look21)


def combine_conj_lookup(look12: np.ndarray, look21: np.ndarray) -> np.ndarray:
    """Pad both tables to square by repeating their last row, then max(look12, look21')."""
    m, n = look12.shape
    look12 = np.vstack([look12, np.repeat(look12[-1:], n - m, axis=0)])
    look21 = np.vstack([look21, np.repeat(look21[-1:], n - m, axis=0)])
    return np.maximum(look12, look21.T)


def cond_fdr(
    lp1: np.ndarray,
    lp2: np.ndarray,
    opts: Options,
    pruneidx: np.ndarray | None = None,
    excludevec: np.ndarray | None = None,
) -> tuple[np.ndarray, tuple, np.ndarray]:
    """Conditional FDR of lp1 given lp2 and the unconditional FDR of lp1 (cond_FDR_amd.m).

    Returns (fdrvec, (looktable, lookcount, lookstd), fdrvec0).
    """
    lp1_tmp, lp2_tmp = np.array(lp1, dtype=float), np.array(lp2, dtype=float)
    if excludevec is not None and np.size(excludevec):
        lp1_tmp[excludevec] = np.nan  # excluded SNPs do not shape the lookup table
        lp2_tmp[excludevec] = np.nan
    if pruneidx is not None and np.size(pruneidx) == 0:
        pruneidx = None
    lookup = lookup_table(lp1_tmp, lp2_tmp, opts, pruneidx)
    looktable = lookup[0]

    hv, t2 = opts.hv, opts.t2breaks
    with np.errstate(invalid="ignore"):
        indvec1 = np.fmin(hv.size, np.fmax(1, 1 + lp1 / (hv[1] - hv[0])))
        indvec2 = np.fmin(opts.t2_nbreaks, np.fmax(1, 1 + lp2 / (t2[1] - t2[0])))
    fdrvec0 = interp2(looktable, indvec1, np.ones_like(indvec1))
    fdrvec = interp2(looktable, indvec1, indvec2)
    fdrvec[~np.isfinite(lp1 + lp2)] = np.nan
    return fdrvec, lookup, fdrvec0
