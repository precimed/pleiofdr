"""CSV result tables (ports of save_to_csv.m, save_fdr.m, save_zscore.m)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .options import Options


@dataclass
class Results:
    fdrmat: np.ndarray
    imat: np.ndarray
    imat2: np.ndarray
    logfdrmat: np.ndarray
    locusnumvec: np.ndarray
    logpmat2: np.ndarray
    logpvec: np.ndarray
    zvec: np.ndarray
    zmat2: np.ndarray
    fdrvec0: np.ndarray
    snpidlist: list[str]
    genenamelist: list[str]
    chrnumvec: np.ndarray
    posvec: np.ndarray
    a1vec: list[str]
    a2vec: list[str]
    flp: np.ndarray
    pruneidx: np.ndarray | None


def _special(x: float) -> str | None:
    if np.isnan(x):
        return "NaN"
    if np.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return None


def fmt_e(x: float) -> str:
    """MATLAB fprintf('%e')."""
    return _special(x) or f"{x:e}"


def fmt_f(x: float) -> str:
    """MATLAB fprintf('%f')."""
    return _special(x) or f"{x:f}"


def fmt_d(x: float) -> str:
    """MATLAB fprintf('%d'): integers as such, non-integers fall back to %e."""
    x = float(x)
    special = _special(x)
    if special:
        return special
    return f"{int(x)}" if x == int(x) else f"{x:e}"


def nanmin(values: np.ndarray) -> float:
    """MATLAB min: ignores NaN, NaN if every value is NaN."""
    finite = values[~np.isnan(values)]
    return float(finite.min()) if finite.size else np.nan


def _stem(outdir: Path, traitname1: str, traitnames: list[str]) -> str:
    return f"{outdir}/{traitname1}{''.join(f'_{t}' for t in traitnames)}"


def _rows(results: Results, pruned: bool) -> np.ndarray:
    selected = np.isfinite(results.locusnumvec)
    if pruned:
        selected &= results.imat2.sum(axis=1) > 0
    return np.flatnonzero(selected)


def _pvals(logfdr: np.ndarray) -> np.ndarray:
    return 10.0 ** (-logfdr)


def save_fdr(
    results: Results,
    traitname1: str,
    traitnames: list[str],
    opts: Options,
    outdir: Path,
    pruned: bool = True,
) -> Path:
    suffix = "loci" if pruned else "all"
    st = opts.stattype
    fname = Path(f"{_stem(outdir, traitname1, traitnames)}_{st}_{opts.fdrthresh:g}_{suffix}.csv")
    header = f"locusnum,snpid,geneid,chrnum,chrpos,pval_{traitname1},fdr_{traitname1}"
    header += "".join(f",{st}_{traitname1}_{t},prune_{traitname1}_{t}" for t in traitnames)
    header += f",min_{st}"
    lines = [header]
    for i in _rows(results, pruned):
        row = [
            fmt_d(results.locusnumvec[i]),
            results.snpidlist[i],
            results.genenamelist[i],
            fmt_d(results.chrnumvec[i]),
            fmt_d(results.posvec[i]),
            fmt_e(10.0 ** -results.logpvec[i]),
            fmt_e(results.fdrvec0[i]),
        ]
        for j in range(len(traitnames)):
            row += [fmt_e(results.fdrmat[i, j]), fmt_d(results.imat2[i, j])]
        row.append(fmt_e(nanmin(_pvals(results.logfdrmat[i, :]))))
        lines.append(",".join(row))
    fname.write_text("\n".join(lines) + "\n")
    return fname


def save_zscore(
    results: Results,
    traitname1: str,
    traitnames: list[str],
    opts: Options,
    outdir: Path,
    pruned: bool = True,
) -> Path:
    suffix = "loci" if pruned else "all"
    st = opts.stattype
    stem = _stem(outdir, traitname1, traitnames)
    fname = Path(f"{stem}_zscore_{st}_{opts.fdrthresh:g}_{suffix}.csv")
    header = f"locusnum,snpid,geneid,chrnum,chrpos,A1,A2,zscore_{traitname1}"
    header += "".join(f",zscore_{t}" for t in traitnames)
    header += "".join(f",{st}_{traitname1}_{t},prune_{traitname1}_{t}" for t in traitnames)
    header += f",min_{st},pval_{traitname1}"
    header += "".join(f",pval_{t}" for t in traitnames)
    lines = [header]
    ntraits = len(traitnames)
    for i in _rows(results, pruned):
        row = [
            fmt_d(results.locusnumvec[i]),
            results.snpidlist[i],
            results.genenamelist[i],
            fmt_d(results.chrnumvec[i]),
            fmt_d(results.posvec[i]),
            results.a1vec[i],
            results.a2vec[i],
            fmt_f(results.zvec[i]),
        ]
        row += [fmt_f(results.zmat2[i, j]) for j in range(ntraits)]
        for j in range(ntraits):
            row += [fmt_e(10.0 ** -results.logfdrmat[i, j]), fmt_d(results.imat2[i, j])]
        row.append(fmt_e(nanmin(_pvals(results.logfdrmat[i, :]))))
        row.append(fmt_e(10.0 ** -results.logpvec[i]))
        row += [fmt_e(10.0 ** -results.logpmat2[i, j]) for j in range(ntraits)]
        lines.append(",".join(row))
    fname.write_text("\n".join(lines) + "\n")
    return fname


def save_to_csv(
    results: Results,
    opts: Options,
    traitname1: str,
    traitnames: list[str],
    outputdir: str | Path,
    prunecsv: bool = True,
) -> list[Path]:
    outdir = Path(outputdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = [save_fdr(results, traitname1, traitnames, opts, outdir, prunecsv)]
    if opts.stattype == "conjfdr":
        written.append(save_zscore(results, traitname1, traitnames, opts, outdir, prunecsv))
    return written
