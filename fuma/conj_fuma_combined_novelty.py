import numpy as np
import pandas as pd
import glob
import os
import sys
from pandas import ExcelWriter

print('This script will take 7 arguments to produce conjFDR_0.05_trait1_vs_trait2_novelty.xlsx')
print('1. combined table of conjFDR_lead & FUMA')
print('2. gwascatalog.txt from respective FUMA out,3. combined table of conjFDR_snp & FUMA 4. Inhouse Novelty db')
print('5. String to search in gwascatalog; i.e. "Depress" 6. TRAIT1 & 7. TRAIT2')


def novelty_check(lead,gwas,snp,db,nov,trait1,trait2):
    nv1=pd.read_csv(lead)
    gw1=pd.read_csv(gwas,sep='\t')
    snp1=pd.read_csv(snp)
    db1=pd.read_csv(db,sep='\t')
    #nov1=nov
    #trait1=sys.argv[6]
    #trait2=sys.argv[7]
    #lead,gwas,snp,db,nov,trait1,trait2=nv1,gw1,snp1,db1,nov1,traits1,trait2
    gw2=gw1[gw1['Trait'].str.contains(fr'{nov}(?!$)')]
    snp1['Novel_in_GWAScatalog'] = snp1.groupby('locusnum')['CAND_SNP'].apply(lambda s: s.isin(gw2['snp']))
    nv1=nv1.merge(snp1.groupby('locusnum')['Novel_in_GWAScatalog'].any(), left_on='locusnum', right_index=True, how='left')
    nv1.Novel_in_GWAScatalog=nv1.Novel_in_GWAScatalog.map({True: 'No', False: 'Yes'})
    nv1['MinBP']=nv1['MinBP'].astype(int)
    nv1['MaxBP']=nv1['MaxBP'].astype(int)
    nv1.CHR=nv1.CHR.astype(int)
    db1=db1[db1.chromosome.str.isnumeric()]
    db1.chromosome=db1.chromosome.astype(int)
    db1['min_bp']=db1['min_bp'].astype(float)
    db1['max_bp']=db1['max_bp'].astype(int)
    a = nv1['MinBP'].to_numpy()
    b = nv1['MaxBP'].to_numpy()
    c1= nv1.CHR.to_numpy()
    c2= db1.chromosome.to_numpy()[:,None]
    # Raise the columns in db1 by a dimension to enable numpy's broadcasting
    c = db1['min_bp'].to_numpy()[:, None]
    d = db1['max_bp'].to_numpy()[:, None]
    result = ((c1 == c2) & (a >= c) & (a <= d)) | ((c1 == c2) & (b >= c) & (b <= d))
    nv1['Novel_in_Database'] = ~result.any(axis=0)
    nv1.Novel_in_Database=nv1.Novel_in_Database.map({True: 'Yes', False: 'No'})
    nv1[f'Novel_in_{trait1}']=np.where(((nv1.Novel_in_GWAScatalog=="Yes") & (nv1.Novel_in_Database=='Yes')),'Yes','No')
    nv1=nv1.drop(['Novel_in_GWAScatalog', 'Novel_in_Database'], axis = 1)
    nv1.to_csv(f'conjFDR_0.05_{trait1}_vs_{trait2}_novelty.csv',index=False)
    
novelty_check(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])

