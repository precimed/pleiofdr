"""Pleiotropy-informed conditional / conjunctional FDR (port of pleioFDR_amd.m).

Reference: Andreassen OA et al. (2013), Am J Hum Genet 92(2):197-209.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .ldmatrix import LDMatrix
from .lookup import cond_fdr
from .options import Options
from .pruning import random_prune_idx


@dataclass
class PleioFDRResult:
    fdrmat: np.ndarray  # (nsnp, ntraits) cond- or conjFDR
    fdrvec0: np.ndarray  # unconditional FDR of trait 1, from unpruned SNPs
    pruneidx: np.ndarray | None
    fdrmat12: np.ndarray  # trait 1 | trait 2
    fdrmat21: np.ndarray  # trait 2 | trait 1 (conjfdr only)
    lookup12: list[tuple]
    lookup21: list[tuple]


def pleio_fdr(
    logpvec1: np.ndarray,
    logpmat2: np.ndarray,
    opts: Options,
    ld: LDMatrix,
    excludevec: np.ndarray,
    pruneidx: np.ndarray | None = None,
) -> PleioFDRResult:
    nsnp, ncondtraits = logpmat2.shape
    fdrmat12 = np.full((nsnp, ncondtraits), np.nan)
    fdrmat21 = np.full((nsnp, ncondtraits), np.nan)
    lookup12: list[tuple] = []
    lookup21: list[tuple] = []
    conj = opts.stattype == "conjfdr"

    if opts.randprune:
        if pruneidx is None or pruneidx.size == 0:
            defvec = ~excludevec & np.isfinite(logpvec1 + logpmat2.sum(axis=1))
            rng = np.random.default_rng(opts.seed)
            # random_prune_idx_amd.m: no post-processing of repeated SNPs
            pruneidx = random_prune_idx(opts.randprune_n, ld, defvec, "default", rng)
    else:
        pruneidx = None

    for i in range(ncondtraits):
        print(f"\n   Trait {i + 1}/{ncondtraits}... ", end="")
        fdrmat12[:, i], look, fdrvec0 = cond_fdr(
            logpvec1, logpmat2[:, i], opts, pruneidx, excludevec
        )
        lookup12.append(look)
        if conj:
            fdrmat21[:, i], look, _ = cond_fdr(logpmat2[:, i], logpvec1, opts, pruneidx, excludevec)
            lookup21.append(look)

    if opts.randprune:
        # the unconditioned fdr is calculated from unpruned SNPs
        _, _, fdrvec0 = cond_fdr(logpvec1, logpmat2[:, 0], opts, None, None)

    fdrmat = np.fmax(fdrmat12, fdrmat21) if conj else fdrmat12
    return PleioFDRResult(fdrmat, fdrvec0, pruneidx, fdrmat12, fdrmat21, lookup12, lookup21)
