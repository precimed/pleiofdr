"""Statistical helpers: binofit, genomic control, Fisher combination, sample overlap."""

from __future__ import annotations

import numpy as np
from scipy import special, stats


def logp_to_absz(logp: np.ndarray) -> np.ndarray:
    """|z| from -log10(p), MATLAB abs(norminv(1/2*10.^-logp))."""
    return np.abs(special.ndtri(0.5 * 10.0 ** (-np.asarray(logp, dtype=float))))


def z_from_logp(logp: np.ndarray) -> np.ndarray:
    """-norminv(0.5 * 10.^-logp) (positive for significant p)."""
    return -special.ndtri(0.5 * 10.0 ** (-np.asarray(logp, dtype=float)))


def logp_from_z(z: np.ndarray) -> np.ndarray:
    """-log10(2*normcdf(-|z|)); p = 0 gives Inf, as in MATLAB."""
    with np.errstate(divide="ignore"):
        return -np.log10(2.0 * special.ndtr(-np.abs(z)))


def binofit(x: np.ndarray, n: np.ndarray, alpha: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Clopper-Pearson estimate and interval, as MATLAB binofit/statbinoci (F-distribution form)."""
    x = np.asarray(x, dtype=float)
    n = np.broadcast_to(np.asarray(n, dtype=float), x.shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        phat = x / n
        nu1, nu2 = 2 * x, 2 * (n - x + 1)
        f = stats.f.ppf(alpha / 2, nu1, nu2)
        lb = (nu1 * f) / (nu2 + nu1 * f)
        nu1, nu2 = 2 * (x + 1), 2 * (n - x)
        f = stats.f.ppf(1 - alpha / 2, nu1, nu2)
        ub = (nu1 * f) / (nu2 + nu1 * f)
    lb[x == 0] = 0.0
    ub[x == n] = 1.0
    return phat, np.column_stack([lb.ravel(), ub.ravel()])


def matlab_prctile(values: np.ndarray, percents: np.ndarray) -> np.ndarray:
    """MATLAB prctile on a vector: NaN ignored, midpoint ('hazen') interpolation."""
    values = values[~np.isnan(values)]
    if values.size == 0:
        return np.full(np.shape(percents), np.nan)
    return np.percentile(values, percents, method="hazen")


def gc_correct_logp(
    logpvec: np.ndarray,
    ivec0: np.ndarray,
    use_standard_gc: bool = False,
    pruneidx: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Genomic control of -log10(p) (GCcorrect_logpvec.m; Devlin & Roeder 1999).

    The in-house estimator takes the median over 50 upper-tail quantiles of z^2 / chi2(1),
    instead of the single median quantile.
    """
    ivec0 = np.asarray(ivec0, dtype=bool)
    if pruneidx is None or np.size(pruneidx) == 0:
        pruneidx = np.ones((logpvec.shape[0], 1), dtype=bool)
    pruneidx = np.asarray(pruneidx, dtype=bool).reshape(logpvec.shape[0], -1)

    fracvec = 1 - np.logspace(np.log10(0.5), -2, 50)
    sig0 = np.full(pruneidx.shape[1], np.nan)
    for i in range(pruneidx.shape[1]):
        zvec0 = logp_to_absz(logpvec[ivec0 & pruneidx[:, i]])
        if use_standard_gc:
            z2 = zvec0[~np.isnan(zvec0)] ** 2
            med = np.median(z2) if z2.size else np.nan
            sig0[i] = np.sqrt(med / stats.chi2.ppf(0.5, 1))
        else:
            with np.errstate(invalid="ignore"):
                prc = matlab_prctile(zvec0**2, 100 * fracvec)
            sig0[i] = np.median(np.sqrt(prc / stats.chi2.ppf(fracvec, 1)))

    zvec = logp_to_absz(logpvec)
    return logp_from_z(zvec / np.median(sig0)), sig0


def fisher_combined(lp1: np.ndarray, lp2: np.ndarray) -> np.ndarray:
    """Fisher combined -log10(p) of trait 1 with each column of lp2 (fisher_comStats.m).

    fisher_comStats.m only GC-corrects when ivec0 is missing, which never happens in the pipeline,
    so no correction is applied here (see MIGRATION_NOTES.md).
    """
    lp2 = lp2.reshape(lp2.shape[0], -1)
    z1sq = special.ndtri(10.0 ** (-lp1) / 2) ** 2
    out = np.empty(lp2.shape)
    for j in range(lp2.shape[1]):
        chi = (z1sq + special.ndtri(10.0 ** (-lp2[:, j]) / 2) ** 2) / 2
        with np.errstate(divide="ignore"):
            out[:, j] = -np.log10(special.gammaincc(1.0, chi))
    return out


def _zscores(logpvec1: np.ndarray, logpmat2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return z_from_logp(logpvec1), z_from_logp(logpmat2)


def check_sample_overlap(
    logpvec1: np.ndarray,
    logpmat2: np.ndarray,
    ivec0: np.ndarray,
    traitname1: str,
    traitnames: list[str],
) -> list[float]:
    print("Testing for sample overlap...")
    z1, z2 = _zscores(logpvec1, logpmat2)
    correlations = []
    for j in range(logpmat2.shape[1]):
        idx = ivec0 & ~np.isnan(z1) & ~np.isnan(z2[:, j])
        c = float(np.corrcoef(z1[idx], z2[idx, j])[0, 1])
        correlations.append(c)
        print(f"\tCorrelation between z scores in {traitname1} and {traitnames[j]} is {c:.3f}")
    print("\tNote that correlation is calculated across intergenic SNPs. A large correlation")
    print("\tmay indicate sample overlap.")
    return correlations


def correct_sample_overlap(
    logpvec1: np.ndarray,
    logpmat2: np.ndarray,
    ivec0: np.ndarray,
    traitname1: str,
    traitnames: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Decorrelate z-scores with C^(-1/2) (Mahalanobis transformation)."""
    if logpmat2.shape[1] != 1:
        raise ValueError("unable to control for sample overlap with more than 2 traits")
    z1, z2 = _zscores(logpvec1, logpmat2[:, 0])
    idx = ivec0 & ~np.isnan(z1) & ~np.isnan(z2)
    c = np.corrcoef(z1[idx], z2[idx])
    print(
        f"Correct correlation between z scores in {traitname1} and {traitnames[0]}; "
        f"before correction the correlation was {c[0, 1]:.3f}"
    )
    evals, evecs = np.linalg.eigh(c)
    c_inv_sqrt = evecs @ np.diag(evals**-0.5) @ evecs.T
    z = c_inv_sqrt @ np.vstack([z1, z2])
    return logp_from_z(z[0]), logp_from_z(z[1])[:, None]
