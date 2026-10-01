"""Independent loci and locus numbering (ports of ind_loci_idx.m and locusnumber.m)."""

from __future__ import annotations

import numpy as np
from numba import njit

from ._caches import NUMBA_CACHE
from .ldmatrix import LDMatrix
from .options import Options
from .pruning import fast_prune


def maf_mask(mafvec: np.ndarray, mafthresh: float) -> np.ndarray:
    """SNPs failing the MAF filter (undefined MAF counts as failing)."""
    if np.isnan(mafthresh):
        return np.zeros(mafvec.shape, dtype=bool)
    return np.isnan(mafvec) | (mafvec <= mafthresh)


def independent_loci(
    fdrmat: np.ndarray, flp: np.ndarray, ld: LDMatrix, mafvec: np.ndarray, opts: Options
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SNPs passing the FDR and p thresholds (imat), LD-pruned lead SNPs (imat2), -log10(FDR)."""
    defvec = np.isfinite(fdrmat + flp)
    with np.errstate(invalid="ignore"):
        imat = (fdrmat <= opts.fdrthresh) & (10.0 ** (-flp) <= opts.pthresh)
    with np.errstate(divide="ignore"):
        logfdrmat = -np.log10(fdrmat)
    logfdrmat[~defvec] = np.nan
    # logfdrmat(mask) = NaN with an nsnp-long mask indexes linearly: first column only
    logfdrmat[maf_mask(mafvec, opts.mafthresh), 0] = np.nan
    logfdrmat_pruned = np.where(imat, logfdrmat, np.nan)

    imat2 = np.zeros(imat.shape, dtype=bool)
    for j in range(fdrmat.shape[1]):
        column = logfdrmat_pruned[:, j]
        tmp = fast_prune(column, ld)
        # only one hit per excluded region survives loci pruning
        for first, last in opts.exclude_from_fit:
            region = column[first : last + 1]
            tmp[first : last + 1] = np.nan
            if np.any(~np.isnan(region)):
                best = first + int(np.nanargmax(region))
                tmp[best] = column[best]
        imat2[:, j] = np.isfinite(tmp)
    return imat, imat2, logfdrmat


@njit(cache=NUMBA_CACHE)
def _locus_walk(ivec: np.ndarray, iivec: np.ndarray, indptr: np.ndarray, indices: np.ndarray):
    n = iivec.size
    locusnumvec = np.full(n, np.nan)
    locusnum = 0
    k = 0  # position in ivec of the first SNP of the next locus
    while k < ivec.size:
        locusnum += 1
        mini = ivec[k]
        maxi = mini
        # extend the locus while SNPs passing threshold are in n-degree LD further right
        while True:
            ii = maxi
            last = -1
            for p in range(indptr[ii], indptr[ii + 1]):
                r = indices[p]
                if iivec[r] and r > last:
                    last = r
            if last <= ii:
                break
            maxi = last
        while k < ivec.size and ivec[k] <= maxi:
            locusnumvec[ivec[k]] = locusnum
            k += 1
    return locusnumvec, locusnum


def locus_number(
    imat: np.ndarray, ld: LDMatrix, mafvec: np.ndarray, opts: Options
) -> tuple[np.ndarray, float]:
    """Label SNPs passing threshold with the number of the LD-connected locus they belong to."""
    iivec = imat.reshape(imat.shape[0], -1).any(axis=1)
    iivec &= ~maf_mask(mafvec, opts.mafthresh)
    ivec = np.flatnonzero(iivec)
    if ivec.size == 0:
        return np.full(iivec.shape, np.nan), np.nan
    locusnumvec, locusnum = _locus_walk(ivec, iivec, ld.indptr, ld.indices)
    return locusnumvec, float(locusnum)
