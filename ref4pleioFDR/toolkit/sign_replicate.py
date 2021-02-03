import pandas as pd
import scipy.io as sio
import numpy as np
from scipy.stats import binom_test

# load reference
ref = pd.read_csv('SUMSTAT/misc/9545380_ref/9545380.ref', sep='\t')
for col in 'CHR BP A1 A2'.split():
    ref.rename(columns={col:'REF_'+col}, inplace=True)
del ref['GP'] 

# add z-scores from matlab files
fnames = '23andMe_MDD_2016 CLOZUK_SCZ_2018_withPGC PGC_MDD_2018_Howard_no23andMe  IHGC_MIG_2016_with23andMe PGC_SCZ_2021_noCLOZUKorPGC2 HUNT_MIG_2021 UKB_MIG_2021 HUNT_MIG_2021_with_UKB_MIG_2021'.split()
pattern = 'SUMSTAT/TMP/mat_9545380/{}.mat'
for fname in fnames:
    mat = sio.loadmat(pattern.format(fname))
    ref[fname] = mat['zvec']
    
# load loci with lead SNPs    
df_mig_scz = pd.read_csv('conj_MIG_SCZ_0.05_results_clump.lead2.csv', sep=',')
df_mig_scz = pd.merge(df_mig_scz, ref, left_on='LEAD_SNP', right_on='SNP', how='left')
for c in '23andMe_MDD_2016 PGC_MDD_2018_Howard_no23andMe'.split(): del df_mig_scz[c]
print(df_mig_scz.shape);

# calc replication status
mig_rep = 'HUNT_MIG_2021_with_UKB_MIG_2021'   # other choices: 'HUNT_MIG_2021' or 'UKB_MIG_2021'
df_mig_scz['SCZ_REP'] = np.sign(df_mig_scz['CLOZUK_SCZ_2018_withPGC'].values) == np.sign(df_mig_scz['PGC_SCZ_2021_noCLOZUKorPGC2'].values)
df_mig_scz['MIG_REP'] = np.sign(df_mig_scz['IHGC_MIG_2016_with23andMe'].values) == np.sign(df_mig_scz[mig_rep].values)

# count how many loci replicated, and compute en-masse replication p-value
yes = np.sum(df_mig_scz['SCZ_REP'].values & df_mig_scz['MIG_REP'].values)
tot = np.sum(np.isfinite(df_mig_scz['PGC_SCZ_2021_noCLOZUKorPGC2'].values) & np.isfinite(df_mig_scz[mig_rep].values))
pval=binom_test(x=tot-yes, n=tot, alternative='less')
print('{} of {} replicate, p-val {}'.format(yes, tot, pval))

# save .csv file
df_mig_scz.to_csv('conj_MIG_SCZ_0.05_results_clump.lead2.repl.csv', sep=',', index=False)
