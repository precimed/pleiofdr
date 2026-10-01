"""Combine pleioFDR clumping output with FUMA annotations and GWAS summary statistics.

Ports of fuma/cond_fuma_combined.R, fuma/conj_fuma_combined_lead.R and
fuma/conj_fuma_combined_snps.R. Column selection is positional, exactly as in the R scripts,
so the inputs must have the layout those scripts were written for.

Usage:
    python -m pleiofdr.fuma.combine {cond,conj-lead,conj-snps} \\
        CLUMP_CSV FUMA_SNPS_TXT SUMSTATS1 SUMSTATS2 TRAIT1 TRAIT2 [--outdir DIR]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

STAT_COLUMNS = ("PVAL", "Z", "OR", "BETA", "SE")

# 1-based column positions kept after the merge, as in the R scripts
LEAD_COLUMNS = [2, 3, 1, 4, 5, 6, 7, 12, 11, 18, 19, 20, 21, 22, 23, 24]
SNPS_COLUMNS = [2, 3, 4, 5, 1, 6, 7, 8, 9, 10, 11, *range(20, 32)]


def read_table(path: str | Path) -> pd.DataFrame:
    """Tab-separated table (gzip detected from the name), integer columns kept integer with NA."""
    return pd.read_csv(
        path,
        sep="\t",
        dtype_backend="numpy_nullable",
        float_precision="round_trip",
        low_memory=False,
    )


def is_true(values: pd.Series) -> pd.Series:
    """is_locus_lead == "True".

    Current data.table::fread parses True/False columns as logical, and the R comparison
    against the string "True" then selects nothing; both spellings are accepted here.
    """
    return values.astype("string").str.lower().eq("true").fillna(False).astype(bool)


def r_merge(left: pd.DataFrame, right: pd.DataFrame, by_x: str, by_y: str) -> pd.DataFrame:
    """R merge(left, right, by.x=, by.y=): key first, then the other left and right columns
    (clashing names get .x/.y suffixes), rows sorted by key."""
    right = right.rename(columns={by_y: by_x})
    lcols = [c for c in left.columns if c != by_x]
    rcols = [c for c in right.columns if c != by_x]
    clash = set(lcols) & set(rcols)
    left = left.rename(columns={c: f"{c}.x" for c in clash})
    right = right.rename(columns={c: f"{c}.y" for c in clash})
    merged = pd.merge(left, right, on=by_x, how="inner", sort=False)
    merged = merged.sort_values(by_x, kind="stable", ignore_index=True)
    order = [by_x] + [c for c in merged.columns if c != by_x]
    return merged[order]


def r_match(keys: pd.Series, table: pd.DataFrame, key: str, value: str) -> pd.Series:
    """table[[value]][match(keys, table[[key]])]: first match, NA when absent."""
    first = table.drop_duplicates(subset=key, keep="first").set_index(key)[value]
    return keys.map(first)


def select_positions(df: pd.DataFrame, positions: list[int]) -> pd.DataFrame:
    return df.iloc[:, [p - 1 for p in positions]].copy()


def add_sumstats(
    df: pd.DataFrame, key: str, sm1: pd.DataFrame, sm2: pd.DataFrame, trait1: str, trait2: str
) -> None:
    """Add <STAT>_in_<trait> columns for statistics present in both summary-statistics files."""
    shared = [s for s in STAT_COLUMNS if s in sm1.columns and s in sm2.columns]
    for sm, trait in ((sm1, trait1), (sm2, trait2)):
        for stat in shared:
            df[f"{stat}_in_{trait}"] = r_match(df[key], sm, "SNP", stat).to_numpy()


def _r_and(a: pd.Series, b: pd.Series) -> pd.Series:
    """R's three-valued &: FALSE & NA is FALSE, TRUE & NA is NA."""
    a, b = a.astype("boolean"), b.astype("boolean")
    return a & b


def _ifelse(cond: pd.Series, yes, no) -> pd.Series:
    """R ifelse(): NA where the condition is NA."""
    cond = cond.astype("boolean")
    yes = yes if isinstance(yes, pd.Series) else pd.Series(yes, index=cond.index)
    no = no if isinstance(no, pd.Series) else pd.Series(no, index=cond.index)
    out = yes.where(cond.fillna(False), no).astype(object)
    out[cond.isna().to_numpy()] = pd.NA
    return out


def add_concordance(
    df: pd.DataFrame, key: str, sm1: pd.DataFrame, sm2: pd.DataFrame, trait1: str, trait2: str
) -> None:
    """Align z-scores to the A1 allele of the table and flag effects with the same sign."""
    df[f"A1_in_{trait1}"] = r_match(df[key], sm1, "SNP", "A1").to_numpy()
    df[f"A1_in_{trait2}"] = r_match(df[key], sm2, "SNP", "A1").to_numpy()
    zrec = {}
    for trait in (trait1, trait2):
        z = pd.to_numeric(df[f"Z_in_{trait}"]).astype("Float64")
        same = df["A1"].astype("string").eq(df[f"A1_in_{trait}"].astype("string"))
        zrec[trait] = pd.to_numeric(_ifelse(same, z, -z)).astype("Float64")
        df[f"Z_recalculated_in_{trait}"] = zrec[trait]
    z1, z2 = zrec[trait1], zrec[trait2]
    df["ConcordEffect"] = _ifelse(
        _r_and(z1 > 0, z2 > 0), "Yes", _ifelse(_r_and(z1 < 0, z2 < 0), "Yes", "No")
    )


def _true_rows(cond: pd.Series) -> pd.Series:
    """subset(): rows where the condition is TRUE (NA dropped)."""
    return cond.astype("boolean").fillna(False).astype(bool)


def combine(
    mode: str,
    clump_csv: str | Path,
    fuma_snps: str | Path,
    sumstats1: str | Path,
    sumstats2: str | Path,
    trait1: str,
    trait2: str,
) -> pd.DataFrame:
    lead = read_table(clump_csv)
    snps = read_table(fuma_snps)
    sm1, sm2 = read_table(sumstats1), read_table(sumstats2)

    if mode == "conj-snps":
        key = "CAND_SNP"
        lead = lead.sort_values(["locusnum", "FDR"], kind="stable", na_position="last")
        lead = lead.drop_duplicates(subset=key, keep="first")
        table = select_positions(r_merge(lead, snps, key, "rsID"), SNPS_COLUMNS)
    else:
        key = "LEAD_SNP"
        lead = lead[is_true(lead["is_locus_lead"])]
        lead = lead.sort_values(["locusnum", "FDR"], kind="stable", na_position="last")
        lead = lead.drop_duplicates(subset="locusnum", keep="first").drop(columns="is_locus_lead")
        table = select_positions(r_merge(lead, snps, key, "rsID"), LEAD_COLUMNS)

    add_sumstats(table, key, sm1, sm2, trait1, trait2)
    if mode != "conj-snps":
        names = list(table.columns)
        names[7], names[8] = "A1", "A2"
        table.columns = names

    fdr = pd.to_numeric(table["FDR"])
    if mode == "cond":
        keep = fdr < 0.01
    elif mode == "conj-lead":
        keep = fdr < 0.05
    else:
        keep = _r_and(fdr < 0.05, pd.to_numeric(table["R2"]) > 0.6)
    table = table[_true_rows(keep)]
    table = table.sort_values("locusnum", kind="stable", na_position="last", ignore_index=True)

    if mode != "cond":
        add_concordance(table, key, sm1, sm2, trait1, trait2)
    return table


def output_name(mode: str, trait1: str, trait2: str) -> str:
    return {
        "cond": f"condFDR_0.01_{trait1}_vs_{trait2}.csv",
        "conj-lead": f"conjFDR_0.05_{trait1}_vs_{trait2}_lead.csv",
        "conj-snps": f"conjFDR_0.05_{trait1}_vs_{trait2}_snps.csv",
    }[mode]


R_DIGITS = 15  # R_print.digits used by write.table ("maximal possible precision")
R_KP_MAX = 22  # powers of ten R keeps in an exact table
R_DEC_MIN_EXPONENT = -308


def _r_scientific(r: float) -> tuple[int, int, bool]:
    """(kpower, nsig, roundingwidens) of |x|, as computed by scientific() in R's format.c.

    R scales |x| to a 15-digit integer (in long double; np.longdouble mirrors the platform),
    drops trailing zeros, and leaves the printed digits to C printf.
    """
    kp = math.floor(math.log10(r)) - R_DIGITS + 1
    r_prec = np.longdouble(r)
    if abs(kp) <= R_KP_MAX:
        r_prec = r_prec / np.longdouble(10.0**kp) if kp >= 0 else r_prec * np.longdouble(10.0**-kp)
    elif kp <= R_DEC_MIN_EXPONENT:
        r_prec = (r_prec * np.longdouble(1e303)) / np.longdouble(10.0 ** (kp + 303))
    else:
        r_prec /= np.longdouble(10.0**kp)
    if r_prec < 10.0 ** (R_DIGITS - 1):
        r_prec *= 10
        kp -= 1
    alpha = float(np.rint(r_prec))  # nearbyint: round half to even
    nsig = R_DIGITS
    for _ in range(R_DIGITS):
        alpha /= 10.0
        if alpha == math.floor(alpha):
            nsig -= 1
        else:
            break
    if nsig == 0:
        nsig, kp = 1, kp + 1
    kpower = kp + R_DIGITS - 1
    rgt = min(max(R_DIGITS - kpower, 0), R_KP_MAX)
    widens = 0 < kpower <= R_KP_MAX and r < 10.0**kpower - 0.5 / 10.0**rgt
    return kpower, nsig, widens


def r_format_number(x: float) -> str:
    """One double as R's write.table prints it: up to 15 significant digits, fixed or scientific
    notation, whichever is narrower (fixed on ties)."""
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0:
        return "0"
    kpower, nsig, widens = _r_scientific(abs(x))
    neg = int(x < 0)
    left = kpower + 1 - int(widens)
    rgt = max(0, nsig - left)
    width_fixed = neg + max(left, 1) + rgt + (rgt != 0)
    exp_digits = 2 if (left > 100 or left <= -99) else 1
    width_sci = neg + (nsig > 1) + (nsig - 1) + 4 + exp_digits
    if width_fixed <= width_sci:
        return f"{x:.{rgt}f}"
    return f"{x:.{nsig - 1}e}"


def _format_cell(value) -> str:
    if value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return "TRUE" if value else "FALSE"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return r_format_number(float(value))
    return str(value)


def write_r_table(df: pd.DataFrame, path: str | Path) -> None:
    """write.table(df, sep='\\t', row.names=FALSE, quote=FALSE)."""
    lines = ["\t".join(map(str, df.columns))]
    for row in df.itertuples(index=False, name=None):
        lines.append("\t".join(_format_cell(v) for v in row))
    Path(path).write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pleiofdr.fuma.combine",
        description="Combine pleioFDR clumping results with FUMA snps.txt and summary statistics",
    )
    parser.add_argument("mode", choices=["cond", "conj-lead", "conj-snps"])
    parser.add_argument("clump_csv", help="clumping output (cond 0.01 / conj 0.05 lead or snps)")
    parser.add_argument("fuma_snps", help="snps.txt from the matching FUMA job")
    parser.add_argument("sumstats1", help="standardised summary statistics for TRAIT1")
    parser.add_argument("sumstats2", help="standardised summary statistics for TRAIT2")
    parser.add_argument("trait1")
    parser.add_argument("trait2")
    parser.add_argument("--outdir", default=".", help="output folder (default: current)")
    args = parser.parse_args(argv)

    table = combine(
        args.mode,
        args.clump_csv,
        args.fuma_snps,
        args.sumstats1,
        args.sumstats2,
        args.trait1,
        args.trait2,
    )
    out = Path(args.outdir) / output_name(args.mode, args.trait1, args.trait2)
    write_r_table(table, out)
    print(f"wrote {out} ({len(table)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
