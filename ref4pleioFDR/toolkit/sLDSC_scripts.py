# Python script to annotate 10M LDSR SNPs into NSS regions (0-1)
import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree

df_annot_nss = pd.read_table(r'Pruefer2014_MM.txt', header=None, names=['CHR', 'FROM', 'TO'])
df_annot_brain = pd.read_csv(r'H:\Dropbox\analysis\2017_02_February_28_cognition_Neanderthal\3.14\hpabraingenes_noheader.txt', sep='\t')

for chri in range(1, 23):
    df = pd.read_csv(r'H:\NORSTORE\MMIL\SUMSTAT\LDSR\LDSR_Annot\1000G_Phase3_baselineLD_ldscores\baselineLD.{0}.annot.gz'.format(chri), delim_whitespace=True)
    df = df[['CHR', 'BP', 'SNP', 'CM']].copy()

    t_nss = df_annot_nss[df_annot_nss['CHR'] == 'chr{}'.format(chri)]
    t_nss = IntervalTree.from_tuples(list(zip(t_nss['FROM'], t_nss['TO'])))
    df['NSS'] = [int(bool(t_nss[p])) for p in df['BP']]

    t_brain = df_annot_brain[df_annot_brain['CHR'] == 'chr{}'.format(chri)]
    t_brain = IntervalTree.from_tuples(list(zip(t_brain['FROM'], t_brain['TO'])))
    df['BRAIN'] = [int(bool(t_brain[p])) for p in df['BP']]

    print(df.shape, df['NSS'].sum(), df['BRAIN'].sum())

    df.to_csv(r'H:\Dropbox\analysis\2018_01_15_NSS\nss.{0}.annot.gz'.format(chri), index=False, sep='\t', compression='gzip')


# Shell script to 
# (1) Calculate LD-weighted annotations
# (2) Execute stratified LD score regression
export LDSCDATA=/mnt/h/NORSTORE/MMIL/SUMSTAT/LDSR/LDSR_Data
export LDSCANNOT=/mnt/h/NORSTORE/MMIL/SUMSTAT/LDSR/LDSR_Annot
export NSSANNOT=/mnt/h/Dropbox/analysis/2018_01_15_NSS
export RESULT=/mnt/h/Dropbox/analysis/2018_01_15_NSS
export TASKS='GIANT_BMI_2015_EUR_lift GIANT_HEIGHT_2014_lift UKB_COLLEGE_2016 SSGAC_EDU_2016 CHARGE_COG_2015 CTG_INTELLIGENCE_2017  UKB_VNR_2016 UKB_RT_2016'
export TASKS='GIANT_BMI_2015_EUR_lift GIANT_HEIGHT_2014_lift CHARGE_COG_2015_lift'
for CHR in {1..22}; do echo "python ldsc.py  --l2 --bfile ${LDSCANNOT}/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHR --ld-wind-cm 1   --print-snps ${LDSCANNOT}/1000G_Phase3_baselineLD_ldscores/list.txt --annot ${NSSANNOT}/nss.$CHR.annot.gz --out ${NSSANNOT}/nss.$CHR &"; done
for TASK in $TASKS; do echo "python ldsc.py --h2 ${LDSCDATA}/${TASK}_noMHC.sumstats.gz --out ${RESULT}/${TASK}.partitioned  --ref-ld-chr ${LDSCANNOT}/1000G_EUR_Phase3_baseline/baseline.,${NSSANNOT}/nss.    --w-ld-chr ${LDSCANNOT}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.   --overlap-annot    --print-coefficients  --frqfile-chr ${LDSCANNOT}/1000G_Phase3_frq/1000G.EUR.QC.    & "; done

# Python script to join all results together
import glob
import os
import re
import pandas as pd
import numpy as np
dir = r'H:\Dropbox\analysis\2018_01_15_NSS\*.partitioned.results'
files = glob.glob(dir)
df_total = None
for fullfile in files:
    file = os.path.split(fullfile)[1]
    df = pd.read_csv(fullfile, delim_whitespace=True)
    df['file'] = file
    df_total = (df if df_total is None else df_total.append(df))
df_total.to_csv(r'H:\Dropbox\analysis\2018_01_15_NSS\partitioned.results.csv', index=False, sep='\t')
