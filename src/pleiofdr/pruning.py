"""LD-based pruning (ports of FastPrune.m, random_prune_idx_amd.m, random_prune_idx_amd_fb.m)."""

from __future__ import annotations

import numpy as np
from numba import njit

from ._caches import NUMBA_CACHE
from .ldmatrix import LDMatrix


@njit(cache=NUMBA_CACHE)
def _prune_mask(order: np.ndarray, indptr: np.ndarray, indices: np.ndarray, n: int) -> np.ndarray:
    prunevec = np.zeros(n, dtype=np.bool_)
    for k in range(order.size):
        s = order[k]
        if not prunevec[s]:
            for p in range(indptr[s], indptr[s + 1]):
                prunevec[indices[p]] = True
            prunevec[s] = False
    return prunevec


def descending_order(values: np.ndarray) -> np.ndarray:
    """Indices of finite |values| in MATLAB sort(abs(x), 'descend') order (stable)."""
    absval = np.abs(values)
    finite = np.flatnonzero(np.isfinite(absval))
    return finite[np.argsort(-absval[finite], kind="stable")]


def fast_prune(logpvec: np.ndarray, ld: LDMatrix) -> np.ndarray:
    """Keep the strongest SNP in each LD block; pruned SNPs are set to NaN.

    Greedy walk from the largest |logp| down: each kept SNP prunes its LD neighbours.
    """
    if logpvec.shape[0] != ld.n:
        raise ValueError("FastPrune: length of logpvec must equal the size of LDmat")
    prunevec = _prune_mask(descending_order(logpvec), ld.indptr, ld.indices, ld.n)
    pruned = np.array(logpvec, dtype=float, copy=True)
    pruned[prunevec] = np.nan
    return pruned


def random_prune_idx(
    niter: int,
    ld: LDMatrix,
    defvec: np.ndarray,
    repeats: str = "maxout",
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Boolean (nsnp, niter) mask of SNPs surviving LD pruning on random priorities."""
    rng = np.random.default_rng() if rng is None else rng
    niter = int(niter)
    prunemask = np.zeros((ld.n, niter), dtype=bool)
    print("\n   Generating random prune indices... ", end="", flush=True)
    for k in range(niter):
        tmp = rng.uniform(0.0, 1.0, ld.n)
        tmp[~defvec] = np.nan  # only SNPs defined in all traits take part
        prunemask[:, k] = np.isfinite(fast_prune(tmp, ld))
        print(f"\r   Generating random prune indices... {k + 1}/{niter}", end="", flush=True)
    print()

    if repeats == "none":
        # at most one iteration per SNP among SNPs selected more than once
        repeated = prunemask.sum(axis=1) > 1
        mask_slice = prunemask[repeated]
        # MATLAB find() enumerates column-major
        cols, rows = np.nonzero(mask_slice.T)
        perm = rng.permutation(rows.size)
        rows, cols = rows[perm], cols[perm]
        _, first = np.unique(rows, return_index=True)
        mask_slice[:] = False
        mask_slice[rows[first], cols[first]] = True
        prunemask[repeated] = mask_slice
    elif repeats == "maxout":
        # cap SNPs present in every iteration at the smallest per-count population, to avoid
        # over-representing singletons
        repeated = prunemask.sum(axis=1) == niter
        repeated_cnt = int(repeated.sum())
        counts = np.bincount(prunemask.sum(axis=1), minlength=niter + 1)[1 : niter + 1]
        mincnt = max(1, int(counts.min()) - 1)
        mask_slice = np.zeros((repeated_cnt, niter), dtype=bool)
        mask_slice[rng.choice(repeated_cnt, mincnt, replace=False), :] = True
        prunemask[repeated] = mask_slice
    elif repeats != "default":
        raise ValueError(f"unknown randprune_repeats option: {repeats}")
    return prunemask
