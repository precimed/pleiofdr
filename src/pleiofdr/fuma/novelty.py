"""Flag conjFDR loci that are novel for trait 1 (port of fuma/conj_fuma_combined_novelty.py).

A locus is novel when none of its candidate SNPs is reported for a matching trait in the
FUMA GWAS catalog and it does not overlap a region of the in-house novelty database.

Usage:
    python -m pleiofdr.fuma.novelty LEAD_TABLE GWASCATALOG_TXT SNPS_TABLE NOVELTY_DB \\
        TRAIT_PATTERN TRAIT1 TRAIT2 [--outdir DIR]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _read_delimited(path: str | Path) -> pd.DataFrame:
    """The combine step writes tab-separated files; comma-separated input is accepted too."""
    with open(path) as f:
        header = f.readline()
    sep = "\t" if "\t" in header else ","
    return pd.read_csv(path, sep=sep, float_precision="round_trip")


def novelty_check(
    lead: str | Path,
    gwas: str | Path,
    snp: str | Path,
    db: str | Path,
    nov: str,
    trait1: str,
) -> pd.DataFrame:
    nv1 = _read_delimited(lead)
    gw1 = pd.read_csv(gwas, sep="\t")
    snp1 = _read_delimited(snp)
    db1 = pd.read_csv(db, sep="\t")

    # catalog traits matching the pattern anywhere except at the very end of the name
    gw2 = gw1[gw1["Trait"].astype(str).str.contains(rf"{nov}(?!$)", regex=True)]
    snp1["Novel_in_GWAScatalog"] = snp1["CAND_SNP"].isin(gw2["snp"])
    in_catalog = snp1.groupby("locusnum")["Novel_in_GWAScatalog"].any()
    nv1 = nv1.merge(in_catalog, left_on="locusnum", right_index=True, how="left")
    nv1["Novel_in_GWAScatalog"] = nv1["Novel_in_GWAScatalog"].map({True: "No", False: "Yes"})

    nv1["MinBP"] = nv1["MinBP"].astype(int)
    nv1["MaxBP"] = nv1["MaxBP"].astype(int)
    nv1["CHR"] = nv1["CHR"].astype(int)
    db1 = db1[db1["chromosome"].astype(str).str.isnumeric()].copy()
    db1["chromosome"] = db1["chromosome"].astype(int)
    db1["min_bp"] = db1["min_bp"].astype(int)
    db1["max_bp"] = db1["max_bp"].astype(int)

    a, b, c1 = nv1["MinBP"].to_numpy(), nv1["MaxBP"].to_numpy(), nv1["CHR"].to_numpy()
    # database regions as rows, loci as columns
    c2 = db1["chromosome"].to_numpy()[:, None]
    c, d = db1["min_bp"].to_numpy()[:, None], db1["max_bp"].to_numpy()[:, None]
    same_chr = c1 == c2
    overlap = same_chr & (((a >= c) & (a <= d)) | ((b >= c) & (b <= d)) | ((a <= c) & (b >= d)))
    nv1["Novel_in_Database"] = np.where(overlap.any(axis=0), "No", "Yes")
    novel = (nv1["Novel_in_GWAScatalog"] == "Yes") & (nv1["Novel_in_Database"] == "Yes")
    nv1[f"Novel_in_{trait1}"] = np.where(novel, "Yes", "No")
    return nv1.drop(columns=["Novel_in_GWAScatalog", "Novel_in_Database"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pleiofdr.fuma.novelty",
        description="Produce conjFDR_0.05_TRAIT1_vs_TRAIT2_novelty.csv",
    )
    parser.add_argument("lead", help="combined conjFDR lead + FUMA table (combine conj-lead)")
    parser.add_argument("gwas", help="gwascatalog.txt from the FUMA job")
    parser.add_argument("snp", help="combined conjFDR snps + FUMA table (combine conj-snps)")
    parser.add_argument("db", help="tab-separated in-house novelty database")
    parser.add_argument("nov", help='text to search for in GWAS catalog traits, e.g. "Depress"')
    parser.add_argument("trait1")
    parser.add_argument("trait2")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args(argv)

    result = novelty_check(args.lead, args.gwas, args.snp, args.db, args.nov, args.trait1)
    out = Path(args.outdir) / f"conjFDR_0.05_{args.trait1}_vs_{args.trait2}_novelty.csv"
    result.to_csv(out, index=False)
    print(f"wrote {out} ({len(result)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
