"""End-to-end self-test on a small synthetic dataset; needs no downloaded data or network.

Usage:
    python -m pleiofdr.selftest [--verbose] [--keep DIR]

Builds a 3,000-SNP reference with block LD and two traits sharing 40 causal variants, runs the
full condFDR and conjFDR analyses and checks that the expected outputs exist and loci are found.
Exits with status 1 on failure.
"""

from __future__ import annotations

import argparse
import contextlib
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp
from scipy import stats

NSNP = 3000
BLOCK = 10  # SNPs per LD block (all pairs within a block in LD)
NCAUSAL = 40


def make_dataset(folder: Path, seed: int = 0) -> None:
    """Write ref.mat, trait1.mat and trait2.mat (MATLAB v5) into folder."""
    rng = np.random.default_rng(seed)
    block_of = np.arange(NSNP) // BLOCK
    same_block = block_of[:, None] == block_of[None, :]
    ldmat = sp.csc_array(same_block)
    chrnumvec = np.where(np.arange(NSNP) < NSNP // 2, 1, 2)
    posvec = 10_000 + 1_000 * np.arange(NSNP)
    scipy.io.savemat(
        folder / "ref.mat",
        {
            "LDmat": ldmat,
            "chrnumvec": chrnumvec[:, None],
            "posvec": posvec[:, None],
            "mafvec": rng.uniform(0.05, 0.5, NSNP)[:, None],
            "is_intergenic": (rng.random(NSNP) < 0.3)[:, None].astype(float),
            "is_ambiguous": np.zeros((NSNP, 1), dtype=bool),
        },
    )
    causal = rng.choice(NSNP // BLOCK, NCAUSAL, replace=False) * BLOCK + BLOCK // 2
    for name in ("trait1", "trait2"):
        z = rng.normal(size=NSNP)
        z[causal] += rng.choice([-1, 1], NCAUSAL) * rng.uniform(4.5, 6.5, NCAUSAL)
        logp = -np.log10(2 * stats.norm.sf(np.abs(z)))
        scipy.io.savemat(folder / f"{name}.mat", {"logpvec": logp[:, None], "zvec": z[:, None]})


def write_config(folder: Path, stattype: str) -> Path:
    config = folder / f"config_{stattype}.txt"
    fdrthresh = 0.05 if stattype == "conjfdr" else 0.01
    config.write_text(
        "\n".join(
            [
                f"reffile={folder / 'ref.mat'}",
                f"traitfolder={folder}",
                "traitfile1=trait1.mat",
                "traitname1=T1",
                "traitfiles={'trait2.mat'}",
                "traitnames={'T2'}",
                f"stattype={stattype}",
                f"fdrthresh={fdrthresh}",
                "randprune_n=3",
                "exclude_chr_pos=[]",
                "onscreen=false",
                f"outputdir={folder / ('out_' + stattype)}",
            ]
        )
        + "\n"
    )
    return config


def check_outputs(outdir: Path, stattype: str) -> int:
    """Raise AssertionError if outputs are missing; return the number of loci found."""
    fdrthresh = "0.05" if stattype == "conjfdr" else "0.01"
    stem = f"T1_T2_{stattype}_{fdrthresh}"
    expected = [
        "result.mat",
        f"{stem}_all.csv",
        f"{stem}_loci.csv",
        f"{stem}_manhattan.png",
        f"{stem}.log",
        "T1_vs_T2_qq.png",
        "T2_vs_T1_lookup.png",
    ]
    if stattype == "conjfdr":
        expected += [f"T1_T2_zscore_{stattype}_{fdrthresh}_loci.csv"]
    missing = [name for name in expected if not (outdir / name).exists()]
    assert not missing, f"{stattype}: missing outputs {missing}"
    nloci = len((outdir / f"{stem}_loci.csv").read_text().splitlines()) - 1
    assert nloci > 0, f"{stattype}: no loci found in the synthetic data"
    fdrmat = scipy.io.loadmat(outdir / "result.mat")["fdrmat"]
    assert fdrmat.shape == (NSNP, 1), f"{stattype}: unexpected fdrmat shape {fdrmat.shape}"
    return nloci


def run(folder: Path, verbose: bool = False) -> dict[str, int]:
    # imported here so that `--help` works quickly
    from .cli import run_from_config

    make_dataset(folder)
    found = {}
    for stattype in ("conjfdr", "condfdr"):
        config = write_config(folder, stattype)
        with contextlib.ExitStack() as stack:
            if not verbose:
                stack.enter_context(
                    contextlib.redirect_stdout(stack.enter_context(open(os.devnull, "w")))
                )
            run_from_config(config, seed=1)
        found[stattype] = check_outputs(folder / f"out_{stattype}", stattype)
    return found


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pleiofdr.selftest", description=__doc__.split("\n\n")[0]
    )
    parser.add_argument("--verbose", action="store_true", help="show the analysis output")
    parser.add_argument("--keep", metavar="DIR", help="write data and results to DIR and keep them")
    args = parser.parse_args(argv)

    try:
        if args.keep:
            folder = Path(args.keep)
            folder.mkdir(parents=True, exist_ok=True)
            found = run(folder, args.verbose)
        else:
            with tempfile.TemporaryDirectory(prefix="pleiofdr-selftest-") as tmp:
                found = run(Path(tmp), args.verbose)
    except Exception as exc:  # report any failure as a failed self-test
        print(f"pleiofdr self-test FAILED: {exc}", file=sys.stderr)
        return 1
    summary = ", ".join(f"{k}: {v} loci" for k, v in found.items())
    print(f"pleiofdr self-test passed ({summary})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
