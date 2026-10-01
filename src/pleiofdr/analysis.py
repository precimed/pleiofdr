"""The pleiotropy analysis pipeline (port of pleiotropy_analysis.m)."""

from __future__ import annotations

import contextlib
import sys
from dataclasses import dataclass, fields
from pathlib import Path
from typing import TextIO

import matplotlib.pyplot as plt
import numpy as np

from . import io
from .fdr import pleio_fdr
from .io import Reference
from .loci import independent_loci, locus_number, maf_mask
from .lookup import cond_fdr
from .options import Options
from .output import Results, fmt_d, save_to_csv
from .plotting import save_figure
from .plotting.lookup import plot_lookup
from .plotting.manhattan import plot_manhattan
from .plotting.qq import plot_enrichment, plot_qq
from .pruning import random_prune_idx
from .stats import (
    check_sample_overlap,
    correct_sample_overlap,
    fisher_combined,
    gc_correct_logp,
    z_from_logp,
)


class _Tee:
    """Write to several streams at once (MATLAB diary)."""

    def __init__(self, *streams: TextIO) -> None:
        self.streams = streams

    def write(self, text: str) -> int:
        for s in self.streams:
            s.write(text)
        return len(text)

    def flush(self) -> None:
        for s in self.streams:
            s.flush()


@contextlib.contextmanager
def diary(path: Path):
    with open(path, "w") as log:
        tee = _Tee(sys.stdout, log)
        with contextlib.redirect_stdout(tee):
            yield


@dataclass
class AnalysisOutput:
    results: Results
    excludevec: np.ndarray
    pruneidx: np.ndarray | None
    fdrmat12: np.ndarray
    fdrmat21: np.ndarray
    lookup12: list
    lookup21: list
    mat_qq: np.ndarray
    mat_qq_inv: np.ndarray
    mat_tdr: np.ndarray
    mat_tdr_inv: np.ndarray
    mat_enrich: np.ndarray
    mat_enrich_inv: np.ndarray
    locusnum: float


def _print_options(opts: Options) -> None:
    print("options =")
    for f in fields(opts):
        print(f"    {f.name:>32}: {getattr(opts, f.name)!r}")


def _names(traitnames: list[str], sep: str = "_", trailing: bool = False) -> str:
    return "".join(f"{t}{sep}" if trailing else f"{sep}{t}" for t in traitnames)


def run_analysis(
    opts: Options,
    ref: Reference,
    traitfile1: str | Path,
    traitfiles: list[str | Path],
    traitname1: str,
    traitnames: list[str],
    refinfo: tuple[list[str], list[str], list[str]] | None = None,
    pruneidx: np.ndarray | None = None,
    make_plots: bool = True,
) -> AnalysisOutput:
    outputdir = Path(opts.outputdir or "test")
    outputdir.mkdir(parents=True, exist_ok=True)
    logname = f"{traitname1}{_names(traitnames)}_{opts.stattype}_{opts.fdrthresh:g}.log"
    with diary(outputdir / logname):
        return _run(
            opts,
            ref,
            traitfile1,
            traitfiles,
            traitname1,
            traitnames,
            refinfo,
            pruneidx,
            outputdir,
            make_plots,
        )


def _run(
    opts: Options,
    ref: Reference,
    traitfile1: str | Path,
    traitfiles: list[str | Path],
    traitname1: str,
    traitnames: list[str],
    refinfo: tuple[list[str], list[str], list[str]] | None,
    pruneidx: np.ndarray | None,
    outputdir: Path,
    make_plots: bool,
) -> AnalysisOutput:
    _print_options(opts)
    ld, nsnp, ivec0, mafvec = ref.ld, ref.nsnp, ref.is_intergenic, ref.mafvec

    # LOAD FILES
    print("Loading GWAS .mat files... ", end="")
    logpvec1, zvec1 = (a[:, 0] for a in io.load_gwas(traitfile1, nsnp, opts.dummy_zscore))
    logpmat2, zmat2 = io.load_gwas(traitfiles, nsnp, opts.dummy_zscore)
    print("done")

    # EXCLUDE SPECIAL SNPS FROM ANALYSIS
    defined = np.isfinite(logpvec1 + logpmat2.sum(axis=1))
    excludevec = np.zeros(nsnp, dtype=bool)
    print(f"{defined.sum()} variants are defined across all traits")
    if not np.isnan(opts.mafthresh):
        maf_excl = maf_mask(mafvec, opts.mafthresh)
        print(
            f"Exclude {maf_excl.sum()} variants due to MAF (minor allele frequency) below "
            f"{opts.mafthresh:.3f} or undefined"
        )
        excludevec |= maf_excl
    for first, last in opts.exclude_from_fit:
        excludevec[first : last + 1] = True  # excluded from the fit, not from discovery
        chrs = np.unique(ref.chrnumvec[first : last + 1])
        chrstr = str(chrs[0]) if chrs.size == 1 else "[" + " ".join(map(str, chrs)) + "]"
        print(
            f"Exclude {last - first + 1} SNPs on chromosome {chrstr} from "
            f"{fmt_d(ref.posvec[first] / 1000)} to {fmt_d(ref.posvec[last] / 1000)} KB"
        )

    excludevec_discovery = np.zeros(nsnp, dtype=bool)
    if opts.exclude_from_fit_and_discovery:
        print("Excluded SNPs will be excluded from both fit and discovery", end="")
        excludevec_discovery = excludevec.copy()
    if opts.exclude_ambiguous_snps:
        excludevec |= ref.is_ambiguous
        excludevec_discovery |= ref.is_ambiguous
        print(
            f"{ref.is_ambiguous.sum()} Ambigous SNPs will be excluded from both fit and discovery"
        )
    if excludevec_discovery.any():
        logpmat2[excludevec_discovery, :] = np.nan
        logpvec1[excludevec_discovery] = np.nan
        zmat2[excludevec_discovery, :] = np.nan
        zvec1[excludevec_discovery] = np.nan
    if not excludevec.any():
        print("No variants were excluded", end="")
    else:
        left = ~excludevec & np.isfinite(logpvec1 + logpmat2.sum(axis=1))
        print(f"{left.sum()} variants left after exclusion described above")

    # RANDOM PRUNING INDICES
    if opts.reset_pruneidx:
        pruneidx = None
    if opts.randprune and pruneidx is None:
        defvec = ~excludevec & np.isfinite(logpvec1 + logpmat2.sum(axis=1))
        rng = np.random.default_rng(opts.seed)
        pruneidx = random_prune_idx(opts.randprune_n, ld, defvec, opts.randprune_repeats, rng)
    if opts.randprune and pruneidx.shape[1] != opts.randprune_n:
        raise ValueError("Invalid pruneidx; try reset_pruneidx=false.")
    if not opts.randprune:
        pruneidx = None

    # GENOMIC CONTROL
    if opts.perform_gc:
        print("\nPerforming genomic correction... ")
        print(
            f"Use {(ivec0 & ~excludevec).sum()} control variants to calculate lambda GC "
            "(genomic correction factor)"
        )
        pruneidx_gc = None
        if opts.randprune_gc and opts.randprune:
            pruneidx_gc = pruneidx
            print("Random pruning is taken into account in lambdaGC calculation. ")
        logpvec1, sig0 = gc_correct_logp(logpvec1, ivec0, opts.use_standard_gc, pruneidx_gc)
        print(f"Genomic Control for {traitname1}: lambda = {np.median(sig0):.3f}")
        zvec1 = z_from_logp(logpvec1) * np.sign(zvec1)
        for i, name in enumerate(traitnames):
            logpmat2[:, i], sig0 = gc_correct_logp(
                logpmat2[:, i], ivec0, opts.use_standard_gc, pruneidx_gc
            )
            print(f"Genomic Control for {name}: lambda = {np.median(sig0):.3f}")
            zmat2[:, i] = z_from_logp(logpmat2[:, i]) * np.sign(zmat2[:, i])
        print("done")
    else:
        print("Skipping genomic correction. ")

    # SAMPLE OVERLAP: the CSV tables and fdrvec0 use the original (non-decorrelated) statistics
    logpvec1_orig, zvec1_orig = logpvec1.copy(), zvec1.copy()
    logpmat2_orig, zmat2_orig = logpmat2.copy(), zmat2.copy()
    if opts.correct_for_sample_overlap:
        logpvec1, logpmat2 = correct_sample_overlap(
            logpvec1, logpmat2, ivec0, traitname1, traitnames
        )
        zvec1 = z_from_logp(logpvec1) * np.sign(zvec1)
        zmat2 = z_from_logp(logpmat2) * np.sign(zmat2)
    else:
        check_sample_overlap(logpvec1, logpmat2, ivec0, traitname1, traitnames)

    # FISHER COMBINED STATISTICS
    print("Calculating Fisher combined stats... ", end="")
    if opts.fishercomb:
        flp = fisher_combined(logpvec1, logpmat2)
    else:
        flp = np.repeat(logpvec1[:, None], logpmat2.shape[1], axis=1)
    flp[np.isnan(logpmat2)] = np.nan
    print("done")

    # PLEIOTROPY ANALYSIS
    print("Running pleiotropy analyses... ", end="")
    pleio = pleio_fdr(logpvec1, logpmat2, opts, ld, excludevec, pruneidx)
    fdrmat = pleio.fdrmat
    # the unconditioned fdr is calculated from unpruned variants
    _, _, fdrvec0 = cond_fdr(logpvec1_orig, logpmat2[:, 0], opts, None, None)
    print("done")
    if pruneidx is not None:
        counts = pruneidx.sum(axis=0)
        std = counts.std(ddof=1) if counts.size > 1 else 0.0
        print(f"mean (std) variants per random pruning iteration = {counts.mean():.2f} ({std:.2f})")

    # PLOT QQ AND ENRICHMENT
    print("PLOT QQ / Enrichment... ", end="")
    args = (logpvec1, logpmat2, traitname1, traitnames, opts, pruneidx, excludevec)
    h_qq, mat_qq = plot_qq(*args)
    h_qq_inv, mat_qq_inv = plot_qq(*args, flip_traits=True)
    h_tdr, mat_tdr = plot_qq(*args, plot_tdr=True)
    h_tdr_inv, mat_tdr_inv = plot_qq(*args, flip_traits=True, plot_tdr=True)
    h_enrich, mat_enrich = plot_enrichment(*args)
    h_enrich_inv, mat_enrich_inv = plot_enrichment(*args, flip_traits=True)

    io.save_result_mat(
        outputdir / "result.mat",
        fdrmat=fdrmat,
        logpvec1=logpvec1,
        logpmat2=logpmat2,
        zvec1=zvec1,
        zmat2=zmat2,
        excludevec=excludevec,
        mafvec=mafvec,
        mat_qq=mat_qq,
        mat_qq_inv=mat_qq_inv,
        mat_enrich=mat_enrich,
        mat_enrich_inv=mat_enrich_inv,
        qq_t1breaks=opts.t1breaks[None, :],
    )

    fwd = f"{outputdir}/{traitname1}_vs_{_names(traitnames, trailing=True)}"
    inv = f"{outputdir}/{_names(traitnames, trailing=True)}vs_{traitname1}_"
    figures = [
        (h_qq, f"{fwd}qq"),
        (h_qq_inv, f"{inv}qq"),
        (h_tdr, f"{fwd}tdr"),
        (h_tdr_inv, f"{inv}tdr"),
        (h_enrich, f"{fwd}enrich"),
        (h_enrich_inv, f"{inv}enrich"),
    ]
    if make_plots:
        for fig, stem in figures:
            for fmt in ("png", "svg"):
                save_figure(fig, fmt, stem)
    print("done")

    # PLOT LOOKUP
    h_lookup = plot_lookup(
        logpvec1, logpmat2, traitname1, traitnames, opts, pruneidx, pleio.lookup12, pleio.lookup21
    )
    if make_plots:
        for fmt in ("png", "svg"):
            save_figure(h_lookup[0] if h_lookup else None, fmt, f"{inv}lookup")

    # SAVE FDR-TABLES AND PLOT MANHATTAN AT A SPECIFIC THRESHOLD
    snpidlist, a1vec, a2vec = refinfo if refinfo is not None else ([""] * nsnp,) * 3
    imat, imat2, logfdrmat = independent_loci(fdrmat, flp, ld, mafvec, opts)
    locusnumvec, locusnum = locus_number(imat, ld, mafvec, opts)
    results = Results(
        fdrmat=fdrmat,
        imat=imat,
        imat2=imat2,
        logfdrmat=logfdrmat,
        locusnumvec=locusnumvec,
        logpmat2=logpmat2_orig,
        logpvec=logpvec1_orig,
        zvec=zvec1_orig,
        zmat2=zmat2_orig,
        fdrvec0=fdrvec0,
        snpidlist=snpidlist,
        genenamelist=[""] * nsnp,
        chrnumvec=ref.chrnumvec,
        posvec=ref.posvec,
        a1vec=a1vec,
        a2vec=a2vec,
        flp=flp,
        pruneidx=pruneidx,
    )

    print("Saving .csv... ", end="")
    save_to_csv(results, opts, traitname1, traitnames, outputdir, False)
    save_to_csv(results, opts, traitname1, traitnames, outputdir, True)
    print("done")

    if opts.manh_plot:
        print("Creating Manhattan plots... ", end="")
        h_manhattan = plot_manhattan(results, traitname1, traitnames, opts)
        formats = ("png",) if opts.stattype == "condfdr" else ("png", "svg")
        stem = f"{outputdir}/{traitname1}{_names(traitnames)}_{opts.stattype}"
        stem += f"_{opts.fdrthresh:g}_manhattan"
        if make_plots:
            for fmt in formats:
                save_figure(h_manhattan, fmt, stem)
        print("done")

    if not opts.onscreen:
        plt.close("all")

    return AnalysisOutput(
        results=results,
        excludevec=excludevec,
        pruneidx=pruneidx,
        fdrmat12=pleio.fdrmat12,
        fdrmat21=pleio.fdrmat21,
        lookup12=pleio.lookup12,
        lookup21=pleio.lookup21,
        mat_qq=mat_qq,
        mat_qq_inv=mat_qq_inv,
        mat_tdr=mat_tdr,
        mat_tdr_inv=mat_tdr_inv,
        mat_enrich=mat_enrich,
        mat_enrich_inv=mat_enrich_inv,
        locusnum=locusnum,
    )
