### These scripts are useful to combine pleioFDR and FUMA results
#### 1. cond_fuma_combined.R
This scripts will take 6 arguments to create condFDR_0.01_TRAIT1_vs_TRAIT2.csv file.  
The arguments are : 1. cond.0.01_clump.lead.csv, 2. respective fuma snps.txt, 3. std-sumstats of trait1, 4. std-sumstats of trait2, 
5. TRAIT1 name 6. TRAIT2 name
##### Example:
Rscript ../cond_fuma_combined.R ../fuma_input/cond.md_crp.lead.csv ../fuma_output_fin1/FUMA_cond_md_crp_lead_job139929/snps.txt ../sumstat-std/PGC_MD_2018_with23andMe_noUKBB.sumstats.gz ../sumstat-std/CHARGE_CRP_2018.sumstats.gz DEP CRP
#### 2. conj_fuma_combined_lead.R
This scripts will take 6 arguments to create conjFDR_0.05_TRAIT1_vs_TRAIT2.csv file. 
The arguments are : 1. conj.0.05_clump.lead.csv, 2. respective fuma snps.txt, 3. std-sumstats of trait1, 4. std-sumstats of trait2, 
5. TRAIT1 name 6. TRAIT2 name
##### Example:
Rscript ../conj_fuma_combined_lead.R ../fuma_input/conj.md_crp.lead.csv ../fuma_output_fin1/FUMA_conj_md_crp_lead_job1399/snps.txt ../sumstat-std/PGC_MD_2018_with23andMe_noUKBB.sumstats.gz ../sumstat-std/CHARGE_CRP_2018.sumstats.gz DEP CRP

#### 3. conj_fuma_combined_snps.R
This script would take 6 arguments to create conjFDR_0.05_TRAIT1_vs_TRAIT2_snps.csv file. 
The arguments are: 1. conj.0.05_clump.snps.csv, 2. respective FUMA snps.tx, 3. std sumstats for TRAIT1
4. std sumstats for TRAIT2, 5. TRAIT1 name, 6. TRAIT2 name
##### Example:
Rscript ../conj_fuma_combined_snps.R ../fuma_input/conj.md_crp.snps.csv ../fuma_output_fin1/FUMA_conj_md_crp_lead_job1399/snps.txt ../sumstat-std/PGC_MD_2018_with23andMe_noUKBB.sumstats.gz ../sumstat-std/CHARGE_CRP_2018.sumstats.gz DEP CRP

#### 4. conj_fuma_combined_novelty.py
This scripts will take 7 arguments to create conjFDR_0.05_TRAIT1_vs_TRAIT2_novelty.csv file.  
The arguments are : 1. conjFDR_0.05_TRAIT1_vs_TRAIT2.csv (created by conj_fuma_combined_lead.R), 2. gwascatalog.txt from respective FUMA out, 3. conjFDR_0.05_TRAIT1_vs_TRAIT2.csv (created by conj_fuma_combined_snps.R) 4. Inhouse Novelty db 5. String to search in gwascatalog; i.e. "Depress" 6. TRAIT1 name & 7. TRAIT2 name.
The Inhouse Novelty is the delimited csv file created based on the corresponding data in "https://drive.google.com/drive/folders/1eagc2z3RdYgIgyudc5u6ru_lN5MScaC_?usp=sharing". For example for depression, copy range A22:H1103 to a new excel file. Save the excel file as txt delimited and then replace .txt with .csv. 
##### Example:
python conj_fuma_combined_novelty.py conjFDR_0.05_DEP_vs_BMI.csv ../fuma_output_fin1/FUMA_conj_md_bmi_lead_job1399/gwascatalog.txt conjFDR_0.05_DEP_vs_BMI_SNPs.csv novelty_db_dep.csv Depress DEP BMI

##### 5. csv_to_excel.ipynb
This notebook is an exploratory code compilation which helped a lot to combine multiple csv files into one excel files with multiple sheets
