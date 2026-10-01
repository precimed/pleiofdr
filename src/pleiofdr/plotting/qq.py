"""Conditional Q-Q, TDR and fold-enrichment plots (plot_qq_amd.m, plot_enrichment_amd.m)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ..lookup import histc
from ..options import Options
from ..stats import binofit
from . import MATLAB_COLORS, subplot_grid, tex


def qq_matrix(
    lp1: np.ndarray,
    lp2: np.ndarray,
    hv: np.ndarray,
    thresholds: np.ndarray,
    pruneidx: np.ndarray | None,
    randprune_n: int,
) -> np.ndarray:
    """Empirical 1 - F(lp1 | lp2 >= threshold) on the grid hv, averaged over prune iterations."""
    qqmat = np.zeros((hv.size, thresholds.size))

    def one_minus_cdf(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
        hc = histc(a[b >= t], hv)
        phat, _ = binofit(np.cumsum(hc), hc.sum())
        return 1 - phat

    if pruneidx is None:
        for j, t in enumerate(thresholds):
            with np.errstate(invalid="ignore"):
                qqmat[:, j] = one_minus_cdf(lp1, lp2, t)
        return qqmat

    cntmat = np.zeros_like(qqmat)
    for k in range(int(randprune_n)):
        keep = pruneidx[:, k]
        a, b = np.where(keep, lp1, np.nan), np.where(keep, lp2, np.nan)
        for j, t in enumerate(thresholds):
            with np.errstate(invalid="ignore"):
                q = one_minus_cdf(a, b, t)
            defvec = np.isfinite(q)
            if not defvec.any():
                continue
            qqmat[defvec, j] += q[defvec]
            cntmat[defvec, j] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return qqmat / cntmat


def _prepare(
    logpvec1: np.ndarray, logpmat2: np.ndarray, excludevec: np.ndarray | None
) -> tuple[np.ndarray, np.ndarray]:
    lp1, lp2 = logpvec1.copy(), logpmat2.copy()
    if excludevec is not None:
        lp1[excludevec] = np.nan
        lp2[excludevec, :] = np.nan
    return lp1, lp2


def _legend_labels(thresholds: np.ndarray, trait2: str) -> list[str]:
    return ["All SNPs" if t == 0 else rf"$p_{{{tex(trait2)}}} < 10^{{-{t:g}}}$" for t in thresholds]


def plot_qq(
    logpvec1: np.ndarray,
    logpmat2: np.ndarray,
    traitname1: str,
    traitnames: list[str],
    opts: Options,
    pruneidx: np.ndarray | None,
    excludevec: np.ndarray | None,
    flip_traits: bool = False,
    plot_tdr: bool = False,
):
    """Q-Q (or TDR) figure; returns (figure, qqmat of the last conditioning trait)."""
    lp1_all, lp2_all = _prepare(logpvec1, logpmat2, excludevec)
    thresholds, hv = opts.qqbreaks, opts.t1breaks
    prune = pruneidx if opts.randprune else None
    nrows, ncols = subplot_grid(len(traitnames))
    fig, axes = plt.subplots(nrows, ncols, figsize=(9 * ncols, 9 * nrows), squeeze=False)
    qqmat = np.zeros((hv.size, thresholds.size))
    for i, trait2_name in enumerate(traitnames):
        lp2 = lp2_all[:, i]
        a, b = (lp2, lp1_all) if flip_traits else (lp1_all, lp2)
        qqmat = qq_matrix(a, b, hv, thresholds, prune, opts.randprune_n)
        trait1, trait2 = (trait2_name, traitname1) if flip_traits else (traitname1, trait2_name)

        ax = axes.flat[i]
        handles = []
        for j in range(thresholds.size):
            color = MATLAB_COLORS[j % len(MATLAB_COLORS)]
            if plot_tdr:
                prev = np.concatenate([[1.0], qqmat[:-1, j]])
                with np.errstate(invalid="ignore", divide="ignore"):
                    tdr = 1 - 10.0 ** (-hv) / prev
                tdr[tdr < 0] = 0
                (h,) = ax.plot(hv, tdr, lw=2, color=color)
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    (h,) = ax.plot(-np.log10(qqmat[:, j]), hv, lw=2, color=color)
            handles.append(h)
        labels = _legend_labels(thresholds, trait2)
        if not plot_tdr:
            (h,) = ax.plot(hv, hv, "k--", lw=1.5)
            handles.append(h)
            labels.append("Expected")
            ax.set_ylim(0, 7.3)
            ax.set_yticks(range(8))
            ax.set_ylabel(rf"Nominal $-\log_{{10}}(p_{{{tex(trait1)}}})$", fontsize=24)
            ax.set_xlabel(rf"Empirical $-\log_{{10}}(q_{{{tex(trait1)}}})$", fontsize=24)
        else:
            ax.set_ylim(0, 1)
            ax.set_xlabel(rf"Nominal $-\log_{{10}}(p_{{{tex(trait1)}}})$", fontsize=24)
            ax.set_ylabel(rf"Conditional $TDR_{{{tex(trait1)}|{tex(trait2)}}}$", fontsize=24)
        ax.set_xlim(0, 7.3)
        ax.set_xticks(range(8))
        ax.tick_params(labelsize=20)
        ax.set_title(f"{trait1} | {trait2}", fontsize=26)
        ax.legend(handles, labels, frameon=False, loc="lower right", fontsize=20)
    for ax in list(axes.flat)[len(traitnames) :]:
        ax.set_visible(False)
    fig.tight_layout()
    return fig, qqmat


def plot_enrichment(
    logpvec1: np.ndarray,
    logpmat2: np.ndarray,
    traitname1: str,
    traitnames: list[str],
    opts: Options,
    pruneidx: np.ndarray | None,
    excludevec: np.ndarray | None,
    flip_traits: bool = False,
):
    """Fold-enrichment figure; returns (figure, enrichment matrix of the last trait)."""
    lp1_all, lp2_all = _prepare(logpvec1, logpmat2, excludevec)
    thresholds, hv = opts.qqbreaks, opts.t1breaks
    prune = pruneidx if opts.randprune else None
    nrows, ncols = subplot_grid(len(traitnames))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12 * ncols, 8 * nrows), squeeze=False)
    enrichmat = np.zeros((hv.size, thresholds.size))
    ymax = 0.0
    for i, trait2_name in enumerate(traitnames):
        lp2 = lp2_all[:, i]
        a, b = (lp2, lp1_all) if flip_traits else (lp1_all, lp2)
        qqmat = qq_matrix(a, b, hv, thresholds, prune, opts.randprune_n)
        trait1, trait2 = (trait2_name, traitname1) if flip_traits else (traitname1, trait2_name)
        with np.errstate(invalid="ignore", divide="ignore"):
            enrichmat = qqmat / qqmat[:, :1]

        ax = axes.flat[i]
        for j in range(thresholds.size):
            ax.plot(hv, enrichmat[:, j], lw=2, color=MATLAB_COLORS[j % len(MATLAB_COLORS)])
        ax.set_xlim(0, 7.3)
        ax.legend(_legend_labels(thresholds, trait2), loc="upper left", fontsize=16)
        ax.set_title(f"{trait1} | {trait2}", fontsize=16)
        # the MATLAB code labels the axes with the traits in this order
        y_pair, x_trait = ((trait1, trait2), trait1) if flip_traits else ((trait2, trait1), trait2)
        ax.set_ylabel(f"Fold Enrichment {y_pair[0]} | {y_pair[1]}", fontsize=24)
        ax.set_xlabel(rf"Nominal $-\log_{{10}}(p_{{{tex(x_trait)}}})$", fontsize=24)
        visible = enrichmat[hv <= 7.3]
        finite = visible[np.isfinite(visible)]
        if finite.size:
            ymax = max(ymax, float(finite.max()))
    for ax in list(axes.flat)[: len(traitnames)]:
        ax.set_ylim(0, ymax * 1.05 if ymax > 0 else 1)
    for ax in list(axes.flat)[len(traitnames) :]:
        ax.set_visible(False)
    fig.tight_layout()
    return fig, enrichmat
