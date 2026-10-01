"""Smoke test of the ref4pleioFDR/toolkit annotation pipeline on a tiny synthetic genome."""

import gzip
import subprocess
import sys

import pandas as pd

from .conftest import ROOT

TOOLKIT = ROOT / "ref4pleioFDR" / "toolkit"


def run(script: str, *args) -> None:
    subprocess.run([sys.executable, str(TOOLKIT / script), *map(str, args)], check=True)


def test_annotation_pipeline(tmp_path):
    # one coding gene on chr1 (+ strand) with two exons, and one noncoding transcript
    known_gene = tmp_path / "knownGene.txt"
    known_gene.write_text(
        "uc1\tchr1\t+\t10000\t20000\t11000\t19000\t2\t10000,18000,\t12000,20000,\tP1\ta1\n"
        "uc2\tchr1\t-\t50000\t52000\t50000\t50000\t1\t50000,\t52000,\tP2\ta2\n"
    )
    run("knownGene2annot.py", known_gene, tmp_path / "annot.txt")

    # bed-like annotation of SNP positions: col 4 = SNP id, col 8 = category (0-based 3 / 7)
    annot = pd.read_csv(tmp_path / "annot.txt", sep="\t", header=None)
    snps = {"rs1": 10500, "rs2": 15000, "rs3": 18500, "rs4": 51000, "rs5": 900000}
    rows = []
    for snp, pos in snps.items():
        for _, a in annot[(annot[5] <= pos) & (annot[6] > pos)].iterrows():
            rows.append(["chr1", pos, pos + 1, snp, a[0], a[1], a[2], a[4]])
    pd.DataFrame(rows).to_csv(tmp_path / "snp_annot.txt", sep="\t", header=False, index=False)
    template = pd.DataFrame([["chr1", p, p + 1, s] for s, p in snps.items()])
    template.to_csv(tmp_path / "template.bed", sep="\t", header=False, index=False)

    run(
        "annot2annomat.py",
        tmp_path / "snp_annot.txt",
        tmp_path / "template.bed",
        tmp_path / "annomat.txt.gz",
    )
    run("uniq_annot.py", tmp_path / "annomat.txt.gz", tmp_path / "uniq.txt.gz")
    uniq = pd.read_csv(tmp_path / "uniq.txt.gz", sep="\t", index_col="SNP")
    assert uniq.loc["rs2", "Intron"] == 1
    assert uniq.loc["rs1", "5UTR"] == 1 and uniq.loc["rs1", "Exon"] == 0  # priority order
    assert uniq.loc["rs4", "NoncodingTranscript"] == 1

    ld_dir = tmp_path / "ld"
    ld_dir.mkdir()
    for chrom in range(1, 23):
        with gzip.open(ld_dir / f"chr{chrom}.schork.r2.ld.gz", "wt") as f:
            f.write(" SNP_A SNP_B R2\n")
            if chrom == 1:  # summed r2 >= 1 is needed to inherit an annotation
                f.write(" rs2 rs5 0.6\n rs5 rs2 0.6\n rs3 rs4 0.5\n")
    run("ld_informed_annot.py", tmp_path / "uniq.txt.gz", ld_dir, tmp_path / "ld.txt.gz")
    ld = pd.read_csv(tmp_path / "ld.txt.gz", sep="\t", index_col="SNP")
    assert ld.loc["rs5", "Intron"] == 1  # inherited through LD with rs2
    assert ld.loc["rs5", "Intergenic"] == 0
    assert ld.loc["rs4", "Exon"] == 0  # r2 0.5 is not enough
