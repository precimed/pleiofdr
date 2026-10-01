### Combining pleioFDR and FUMA results

These tools combine pleioFDR clumping output with FUMA annotations and GWAS summary statistics.
They are part of the `pleiofdr` Python package (`pleiofdr.fuma`) and replace the earlier R scripts
(`cond_fuma_combined.R`, `conj_fuma_combined_lead.R`, `conj_fuma_combined_snps.R`) and
`conj_fuma_combined_novelty.py`, producing the same tables.

As in the R scripts, columns are selected by position, so the inputs must have the usual layout:
the clumping `.lead.csv` / `.snps.csv` files, FUMA `snps.txt`, and standardised, tab-separated
summary statistics with `SNP`, `A1` and (some of) `PVAL`, `Z`, `OR`, `BETA`, `SE` columns.
Summary statistics may be gzipped. All outputs are written to the current folder (or `--outdir`).

#### 1. Conditional FDR lead SNPs → `condFDR_0.01_TRAIT1_vs_TRAIT2.csv`
Arguments: 1. cond 0.01 clump `lead.csv`, 2. matching FUMA `snps.txt`, 3. summary statistics of
TRAIT1, 4. summary statistics of TRAIT2, 5. TRAIT1 name, 6. TRAIT2 name.
```
python -m pleiofdr.fuma.combine cond ../fuma_input/cond.md_crp.lead.csv ../fuma_output_fin1/FUMA_cond_md_crp_lead_job139929/snps.txt ../sumstat-std/PGC_MD_2018_with23andMe_noUKBB.sumstats.gz ../sumstat-std/CHARGE_CRP_2018.sumstats.gz DEP CRP
```

#### 2. Conjunctional FDR lead SNPs → `conjFDR_0.05_TRAIT1_vs_TRAIT2_lead.csv`
Same arguments, with the conj 0.05 clump `lead.csv`. Adds allele-aligned z-scores and a
`ConcordEffect` column (same direction of effect in both traits).
```
python -m pleiofdr.fuma.combine conj-lead ../fuma_input/conj.md_crp.lead.csv ../fuma_output_fin1/FUMA_conj_md_crp_lead_job1399/snps.txt ../sumstat-std/PGC_MD_2018_with23andMe_noUKBB.sumstats.gz ../sumstat-std/CHARGE_CRP_2018.sumstats.gz DEP CRP
```

#### 3. Conjunctional FDR candidate SNPs → `conjFDR_0.05_TRAIT1_vs_TRAIT2_snps.csv`
Same arguments, with the conj 0.05 clump `snps.csv`; keeps candidate SNPs with FDR < 0.05 and R2 > 0.6.
```
python -m pleiofdr.fuma.combine conj-snps ../fuma_input/conj.md_crp.snps.csv ../fuma_output_fin1/FUMA_conj_md_crp_lead_job1399/snps.txt ../sumstat-std/PGC_MD_2018_with23andMe_noUKBB.sumstats.gz ../sumstat-std/CHARGE_CRP_2018.sumstats.gz DEP CRP
```

#### 4. Novelty → `conjFDR_0.05_TRAIT1_vs_TRAIT2_novelty.csv`
Arguments: 1. the table from step 2, 2. `gwascatalog.txt` from the FUMA job, 3. the table from
step 3, 4. the in-house novelty database, 5. text to search for in GWAS catalog traits
(e.g. "Depress"), 6. TRAIT1 name, 7. TRAIT2 name. Tab- and comma-separated tables are both accepted.

The in-house novelty database is a tab-separated file built from the corresponding data in
"https://drive.google.com/drive/folders/1eagc2z3RdYgIgyudc5u6ru_lN5MScaC_?usp=sharing". For
example, for depression, copy range A22:H1103 to a new Excel file and save it as tab-delimited text.
```
python -m pleiofdr.fuma.novelty conjFDR_0.05_DEP_vs_BMI_lead.csv ../fuma_output_fin1/FUMA_conj_md_bmi_lead_job1399/gwascatalog.txt conjFDR_0.05_DEP_vs_BMI_snps.csv novelty_db_dep.csv Depress DEP BMI
```

#### 5. csv_to_excel.ipynb
An exploratory notebook that combines several csv files into one Excel workbook with one sheet per
file (needs `openpyxl`).
