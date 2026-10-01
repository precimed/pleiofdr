#!/bin/bash
# Regenerate the R reference outputs for tests/test_fuma.py from the original R scripts
# (git show d814b8a:fuma/*.R, before they were replaced by pleiofdr.fuma).
# Needs Rscript with the data.table package. Run from this folder.
set -euo pipefail
here=$(pwd)
tmp=$(mktemp -d)
for s in cond_fuma_combined conj_fuma_combined_lead conj_fuma_combined_snps; do
    git show "d814b8a:fuma/$s.R" > "$tmp/$s.R"
done
cd "$tmp"
Rscript cond_fuma_combined.R "$here/lead.csv" "$here/fuma_snps.txt" "$here/sumstats1.txt.gz" "$here/sumstats2.txt" T1 T2
Rscript conj_fuma_combined_lead.R "$here/lead.csv" "$here/fuma_snps.txt" "$here/sumstats1.txt.gz" "$here/sumstats2.txt" T1 T2
Rscript conj_fuma_combined_snps.R "$here/snps.csv" "$here/fuma_snps.txt" "$here/sumstats1.txt.gz" "$here/sumstats2.txt" T1 T2
for f in *.csv; do cp "$f" "$here/expected_$f"; done
