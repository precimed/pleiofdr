"""Command-line entry point (port of runme.m)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from . import __version__, io
from .analysis import run_analysis
from .config import TextConfig
from .options import Options
from .plotting import use_backend


def options_from_config(cfg: TextConfig) -> Options:
    return Options(
        randprune=cfg.get_bool("randprune"),
        randprune_gc=cfg.get_bool("randprune_gc"),
        reset_pruneidx=cfg.get_bool("reset_pruneidx"),
        randprune_n=int(cfg.get_num("randprune_n")),
        randprune_file=cfg.get_str("randprune_file"),
        randprune_repeats=cfg.get_str("randprune_repeats"),
        stattype=cfg.get_str("stattype"),
        fdrthresh=cfg.get_num("fdrthresh"),
        pthresh=cfg.get_num("pthresh"),
        onscreen=cfg.get_bool("onscreen"),
        outputdir=cfg.get_str("outputdir"),
        manh_fontsize_genenames=cfg.get_num("manh_fontsize_genenames"),
        manh_plot=cfg.get_bool("manh_plot"),
        manh_legend=cfg.get_str("manh_legend"),
        manh_yspace=cfg.get_num("manh_yspace"),
        manh_ymargin=cfg.get_num("manh_ymargin"),
        manh_colorlist=0.8 * cfg.get_mat("manh_colorlist"),
        exclude_from_fit_and_discovery=cfg.get_bool("exclude_from_discovery"),
        use_standard_gc=cfg.get_bool("use_standard_gc"),
        perform_gc=cfg.get_bool("perform_gc"),
        exclude_ambiguous_snps=cfg.get_bool("exclude_ambiguous_snps"),
        dummy_zscore=cfg.get_bool("dummy_zscore"),
        mafthresh=cfg.get_num("mafthresh"),
    )


def exclusion_regions(
    exclude_chr_pos: np.ndarray, chrnumvec: np.ndarray, posvec: np.ndarray
) -> list[tuple[int, int]]:
    """Map [CHR BP_from BP_to] rows to 0-based inclusive (first, last) SNP indices."""
    regions = []
    for chrom, bp_from, bp_to in exclude_chr_pos.reshape(-1, 3) if exclude_chr_pos.size else []:
        first = np.flatnonzero((chrnumvec == chrom) & (posvec >= bp_from))
        last = np.flatnonzero((chrnumvec == chrom) & (posvec <= bp_to))
        if first.size == 0 or last.size == 0:
            raise ValueError(f"exclusion region [{chrom:g} {bp_from:g} {bp_to:g}] matches no SNPs")
        regions.append((int(first[0]), int(last[-1])))
    return regions


def run_from_config(config: str | Path, seed: int | None = None, make_plots: bool = True):
    """Run the analysis described by a config file; returns the AnalysisOutput."""
    cfg = TextConfig()
    if Path(config).is_file():
        cfg.load_file(config)
    else:
        print(f'config file "{config}" not found; using defaults', file=sys.stderr)
    cfg.print()

    traitfolder = cfg.get_str("traitfolder")
    traitfile1 = cfg.get_str("traitfile1")
    traitname1 = cfg.get_str("traitname1")
    traitfiles = cfg.get_cell("traitfiles")
    traitnames = cfg.get_cell("traitnames")
    if traitfolder:
        traitfile1 = str(Path(traitfolder) / traitfile1)
        traitfiles = [str(Path(traitfolder) / f) for f in traitfiles]
    if len(traitfiles) != len(traitnames):
        raise ValueError("traitfiles and traitnames must have the same length")

    opts = options_from_config(cfg)
    opts.seed = seed
    use_backend(opts.onscreen)

    reffile = cfg.get_str("reffile")
    print(f'Loading reference file ("{reffile}") ...', end="", flush=True)
    ref = io.load_reference(reffile)
    print("done")

    refinfo = None
    if cfg.get_str("refinfo"):
        path = cfg.get_str("refinfo")
        print(f'Loading additional reference information file ("{path}") ...', end="")
        refinfo = io.load_refinfo(path, ref.nsnp)
        print("done")

    pruneidx = None
    if opts.randprune_file and Path(opts.randprune_file).exists():
        print(f"Loading prune set {opts.randprune_file}... ", end="")
        pruneidx = io.load_pruneidx(opts.randprune_file)
        opts.randprune_n = pruneidx.shape[1]
        opts.reset_pruneidx = False
        print("done")

    opts.exclude_from_fit = exclusion_regions(
        cfg.get_mat("exclude_chr_pos"), ref.chrnumvec, ref.posvec
    )
    opts.validate()

    output = run_analysis(
        opts, ref, traitfile1, traitfiles, traitname1, traitnames, refinfo, pruneidx, make_plots
    )
    if opts.onscreen:
        plt.show()
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pleiofdr",
        description="Pleiotropy-informed conditional and conjunctional false discovery rate",
    )
    parser.add_argument("--config", default="config.txt", help="config file (default: config.txt)")
    parser.add_argument("--seed", type=int, default=None, help="seed for random pruning")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args(argv)
    run_from_config(args.config, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
