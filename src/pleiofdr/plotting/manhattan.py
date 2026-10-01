"""Manhattan plot of -log10(cond/conjFDR) (plot_Manhattan.m).

plot_Manhattan.m hard-codes showgenes = false, so gene labels are not drawn here either.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ..options import Options
from ..output import Results
from . import filter_points

_LEGEND_LOC = {
    "northeast": "upper right",
    "northwest": "upper left",
    "southeast": "lower right",
    "southwest": "lower left",
    "north": "upper center",
    "south": "lower center",
    "east": "center right",
    "west": "center left",
    "best": "best",
}


def plot_manhattan(results: Results, traitname1: str, traitnames: list[str], opts: Options):
    fdrmat, imat, imat2 = results.fdrmat, results.imat, results.imat2
    chrnumvec = results.chrnumvec
    ntraits = fdrmat.shape[1]
    colorlist = np.array(opts.manh_colorlist, dtype=float, copy=True)
    with np.errstate(divide="ignore"):
        logfdrmat = -np.log10(fdrmat)

    if opts.stattype == "condfdr":
        # unconditioned FDR as an extra column, drawn in gray and listed first in the legend
        legends = [f"{traitname1} | {t}" for t in traitnames] + [traitname1]
        fdrvec0 = results.fdrvec0.copy()
        fdrvec0[np.isnan(fdrmat).any(axis=1)] = np.nan
        with np.errstate(divide="ignore"):
            logfdrmat = np.column_stack([logfdrmat, -np.log10(fdrvec0)])
        legendsorder = [ntraits, *range(ntraits)]
        if colorlist.shape[0] < ntraits + 1:
            colorlist = np.vstack([colorlist, np.zeros((ntraits + 1 - colorlist.shape[0], 3))])
        colorlist[ntraits] = [0.5, 0.5, 0.5]
        ylabel = r"$-\log_{10}(condFDR)$"
    else:
        legends = [f"{traitname1} & {t}" for t in traitnames]
        fdrvec0 = results.fdrvec0
        legendsorder = list(range(ntraits))
        ylabel = r"$-\log_{10}(conjFDR)$"

    nsnp = fdrmat.shape[0]
    fig, ax = plt.subplots(figsize=(19.2, 7.2))
    finite = logfdrmat[np.isfinite(logfdrmat)]
    ax.set_xlim(-100000, nsnp)
    ax.set_ylim(-0.1, (finite.max() if finite.size else 0) + 1)

    chromosomes = np.unique(chrnumvec)
    if chromosomes.size > 1:
        ticks, labels = [], []
        for k in chromosomes:
            ind = np.flatnonzero(chrnumvec == k)
            gray = 1 - 0.1 * (1 - k % 2)
            ax.fill_between([ind[0], ind[-1]], 0, 300, color=[gray] * 3, lw=0, zorder=0)
            ticks.append(np.median(ind))
            labels.append(str(k))
        ax.set_xticks(ticks, labels, fontsize=13)
    else:
        ax.set_title(f"Chromosome {chromosomes[0]}", fontsize=13)
        ax.set_xticklabels([])

    handles = []
    for i in legendsorder:
        (h,) = ax.plot(-1, -1, "o", mfc=colorlist[i], mec=colorlist[i], ms=8, ls="none")
        handles.append(h)

    ax.plot([-100000, nsnp - 1], [-np.log10(opts.fdrthresh)] * 2, ":", color="k", lw=1.2)

    indvec = np.arange(nsnp, dtype=float)
    for phase in (1, 2, 3):  # (1) above threshold, (2) significant, (3) loci
        for i in range(logfdrmat.shape[1]):
            if i < ntraits:
                ivec, ivec2 = imat[:, i], imat2[:, i]
            else:
                with np.errstate(invalid="ignore"):
                    ivec = fdrvec0 <= opts.fdrthresh
                ivec2 = np.zeros_like(ivec)
            face = colorlist[i]
            if phase == 1:
                sel, edge, size, lw = ~ivec, face, 4, 1
            elif phase == 2:
                sel, edge, size, lw = ivec & ~ivec2, face, 8, 1
            else:
                sel, edge, size, lw = ivec2, "k", 10, 2
            x, y = filter_points(indvec[sel], logfdrmat[sel, i], opts.manh_image_size)
            ax.plot(x, y, "o", mfc=face, mec=edge, mew=lw, ms=size, ls="none")

    if chromosomes.size > 1:
        loc = _LEGEND_LOC.get(opts.manh_legend.lower(), "upper right")
        ax.legend(handles, [legends[i] for i in legendsorder], frameon=False, loc=loc, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xlabel("Chromosome", fontsize=18)
    ax.tick_params(direction="out", width=1.2)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    fig.tight_layout()
    return fig
