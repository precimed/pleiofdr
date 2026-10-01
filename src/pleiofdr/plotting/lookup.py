"""FDR lookup-table heatmaps (plot_lookup.m)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ..lookup import combine_conj_lookup, conj_lookup_table, lookup_table
from ..options import Options
from . import tex


def plot_lookup(
    logpvec1: np.ndarray,
    logpmat2: np.ndarray,
    traitname1: str,
    traitnames: list[str],
    opts: Options,
    pruneidx: np.ndarray | None,
    lookup12: list[tuple] | None = None,
    lookup21: list[tuple] | None = None,
) -> list:
    """One heatmap per conditioning trait; returns the figures."""
    hv = opts.hv
    t2 = opts.t2breaks
    prune = pruneidx if opts.randprune else None
    figures = []
    for i, trait2 in enumerate(traitnames):
        lp2 = logpmat2[:, i]
        if opts.stattype == "condfdr":
            lookup = lookup12[i][0] if lookup12 else lookup_table(logpvec1, lp2, opts, prune)[0]
        else:
            if lookup12 and lookup21:
                lookup = combine_conj_lookup(lookup12[i][0], lookup21[i][0])
            else:
                lookup = conj_lookup_table(logpvec1, lp2, opts, prune)

        fig, ax = plt.subplots(figsize=(10, 8))
        dx, dy = hv[1] - hv[0], t2[1] - t2[0]
        keep = hv < 7.3
        t1, t2n = tex(traitname1), tex(trait2)
        if opts.stattype == "condfdr":
            extent = (hv[0] - dx / 2, hv[keep][-1] + dx / 2, t2[0] - dy / 2, t2[-1] + dy / 2)
            image = ax.imshow(
                lookup[:, keep],
                origin="lower",
                vmin=0,
                vmax=1,
                cmap="hot_r",
                aspect="auto",
                extent=extent,
                interpolation="nearest",
            )
            ax.set_title(rf"Conditional $FDR_{{{t1}|{t2n}}}$", fontsize=22)
            ax.set_xticks(np.arange(np.ceil(opts.t1_low), hv[keep][-1] + dx / 2, 1))
            ax.set_yticks(np.arange(np.ceil(opts.t2_low), opts.t2_up + dy / 2, 1))
        else:
            sel = lookup[np.ix_(keep, keep)].T
            span = hv[keep][-1] + dx / 2
            image = ax.imshow(
                sel,
                origin="lower",
                vmin=0,
                vmax=1,
                cmap="hot_r",
                extent=(hv[0] - dx / 2, span, hv[0] - dx / 2, span),
                interpolation="nearest",
            )
            ax.set_title(rf"Conjunctional $FDR_{{{t1}&{t2n}}}$", fontsize=22)
            ax.set_xticks(range(1, 8))
            ax.set_yticks(range(1, 8))
        ax.set_xlabel(rf"$-\log_{{10}}(p_{{{t1}}})$", fontsize=22)
        ax.set_ylabel(rf"$-\log_{{10}}(p_{{{t2n}}})$", fontsize=22)
        ax.tick_params(labelsize=22)
        fig.colorbar(image, ax=ax)
        fig.tight_layout()
        figures.append(fig)
    return figures
