#====================================================================================================================================
#
#====Description====
#
#Files for https://github.com/precimed/pleiofdr analysis.
#See https://github.com/precimed/pleiofdr/blob/master/README.md for further details.
#
#Description of the files:
# * 262M May  6 21:45 9545380.ref                        - SNPs used in pleioFDR analysis (SNPs template)
# * 437M May  6 20:37 CTG_COG_2018.mat                   - Summary statistics for Intelligence GWAS [1]
# * 583M May  6 20:21 SSGAC_EDU_2016.mat                 - Summary statistics for Educational attainment GWAS [2]
# * 5.3K May  6 21:51 about.txt                          - this file
# * 2.3G May  6 21:45 ref9545380_1kgPhase3eur_LDr2p1.mat - reference file (LD matrix, etc) in matlab format for pleiofdr analysis
# * 1.4G May  6 21:44 ref9545380_bfile.tar.gz            - 1000 Genomes Phase3 genotypes [3]
#                                                         (converted to plink format and aligned to pleiofdr SNPs template)
# *  21M May  7 20:22 pleioFDR_demo_data.tar.gz          - small package containing demo data (costrained to chr21)  

#[1] Savage JE at el., Genome-wiide association meta-analysis (N=269,867) identifies new genetic and functional links to intelligence. Nature Genetics, 2018 Jul;50(7):912-919
#    https://ctg.cncr.nl/software/summary_statistics , SavageJansen_IntMeta_sumstats.zip

#[2] Okbay et al. (2016), Genome-wide association study identifies 74 loci associated with educational attainment. Nature, 533, 539-542. doi: 10.1038/nature17671
#    https://www.thessgac.org/data , EduYears_Main.txt  (Discovery and replication cohorts except 23andMe)

#[3] A global reference for human genetic variation, The 1000 Genomes Project Consortium, Nature 526, 68-74 (01 October 2015) doi:10.1038/nature15393 
#    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

#The description below give some details on how the reference was generated.

#CTG_COG_2018.mat and SSGAC_EDU_2016.mat were derived from the original summary stats using https://github.com/precimed/python_convert

#====================================================================================================================================
#
#====Summary====
#
# 1.input file: 
# * genotype file from 1000G
# * annotation files from UCSC and Ensembl
# 
# 2. Output file:
# Here are the files you need to use in analysis.
# * a mat file (e.g. ref8452254_1kgPhase3eas_LDr2p1.mat). Data in the ref file includes:
# 	 - chrnumvec
# 	 - id1
# 	 - id2
# 	 - is_ambiguous
# 	 - is_intergenic
# 	 - mafvec 
#  	 - posvec
#	 - LDmat
# * a text file (e.g. 8452254.ref).
#	Header: CHR SNP GP BP A1 A2
# * plink file * chr1-22
#	Those are the ref data for clumping.
#	 - $workdir/dataprocess/chr*bim
#	 - $workdir/dataprocess/chr*bed
#	 - $workdir/dataprocess/chr*fam
#
# 3. requirement
# * Memory requirement
#	Step 2.9 & step4 is memory consuming, make sure your computer or server has memory >128GB.
# * Plink
#	We use plink1.9 here.
# * Python
#	version >3 for step2
#	version >2.7 & <3 for other steps
# * bedtools2
#	install it in $workdir/toolkit
#	git clone https://github.com/arq5x/bedtools2.git
#	make
#
#====================================================================================================================================
#
#====Detailed Scripts====
#
# 0.Define workdir, with three subfolders - rawdata, dataprocess, and toolkit
# * raw data: to save the raw genotype and annotation files.
# * data process: to save the temporary files.
# * tookit: contains the python script and software you need to use.
workdir = /home/ref4ccFDR
mkdir $workdir/rawdata $workdir/dataprocess

# 1. Genotype Preparation 
# 1.1 Download data
# Download data from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ and save it to $workdir/rawdata .
#
# 1.2 Sample selection
# Get sample IDs based on the interested population. (set family ID = individual ID)
# Details for the sample ID and related population can be found in: 
# Here we use east Asian population as an example:
awk 'BEGIN{OFS="\t"; FS="\t"} {if ($7 ~ "CHB|JPT|CHS|CDX|KHV") print($2,$2);}' $workdir/rawdata/integrated_call_samples_v2.20130502.ALL.ped > $workdir/dataprocess/samples.eas
# Extracted 673 samples. After applying additional filters in "plink --vcf" step below, only 504 samples remained.

# Infer sex based on genotypes:
tail -n+2 $workdir/rawdata/integrated_call_samples_v2.20130502.ALL.ped | awk 'BEGIN{OFS="\t"; FS="\t";} {print($2,$2,$5)}' > $workdir/dataprocess/samples.sex

# 1.2. Genotypes
# [1]. Duplicates:
seq 22 | parallel -j 11 "zgrep -v '^#' $workdir/rawdata/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | cut -f 3 | sort | uniq -d > $workdir/dataprocess/chr{}.dups"

# [2]. Convert format(vcf to plink) & set sex based on the genotype & quality control: 
seq 22 | parallel -j 11 "plink --vcf $workdir/rawdata/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --biallelic-only strict --out  $workdir/dataprocess/chr{} --make-bed --mind 0.1 --geno 0.1 --hwe 1.E-20 midp --maf 0.01 --keep $workdir/dataprocess/samples.eas --exclude $workdir/dataprocess/chr{}.dups --update-sex $workdir/dataprocess/samples.sex" &

# [3]. r2:
# Calculate r2:
# Here r2 coefficients are calculated based on genotypes {0,1,2}, to have r2 values based on maximum likelihood haplotypes 'dprime' option should be included
for ((i=1;i<23;i++)); do plink --bfile $workdir/dataprocess/chr${i} --r2 gz yes-really --ld-window 1000000 --ld-window-kb 20000 --ld-window-r2 0.05 --out $workdir/dataprocess/chr${i}.r2; done
# For more details about genotype and haplotype based r2 estimations refer to plink's manual (see 'dprime' option of --r/--r2 flags):
# https://www.cog-genomics.org/plink/1.9/ld#r

# Convert format: plink to space-delimited format (write to tmp dir, than overwrite files in the root 1000genomes dir):
mkdir $workdir/dataprocess/tmp
seq 22 | parallel -j 11 "zcat $workdir/dataprocess/chr{}.r2.ld.gz | sed 's/^[ \t]*//' | tr -s ' ' |  gzip -c > $workdir/dataprocess/tmp/chr{}.r2.ld.gz"

# [4]. MAF:
seq 22 | parallel -j 11 "plink --bfile $workdir/dataprocess/chr{} --freq --out $workdir/dataprocess/chr{}.maf"

# [5]. "identical":
mkdir $workdir/dataprocess/annot
mkdir $workdir/dataprocess/annot/identical
for ((i=1;i<23;i++)); do awk 'BEGIN{OFS="\t"; print("SNP","IDENTICAL");} {print($2,1);}' $workdir/dataprocess/chr${i}.bim > $workdir/dataprocess/annot/identical/chr${i}.annot; done

# [6]. Gather all variants:
echo -e "CHR\tSNP\tBP\tA1\tA2" > $workdir/dataprocess/template9524k.txt && cut -f1,2,4,5,6 $workdir/dataprocess/chr*.bim >> $workdir/dataprocess/template9524k.txt && gzip $workdir/dataprocess/template9524k.txt
less $workdir/dataprocess/template9524k.txt.gz|awk 'NR>1'|awk '{print($1"\t"$2"\t0\t"$3"\t"$4"\t"$5);}'|sort --key=1 -n > $workdir/dataprocess/all_chr.bim &

# 2.Annotation Preparation
# 2.1 Gene info
# Data source: KnownGene from UCSC
wget -O $workdir/rawdata/knownGene_ucsc_hg19.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz

# Extract infoList
python2 $workdir/toolkit/knownGene2annot.py $workdir/rawdata/knownGene_ucsc_hg19.txt.gz $workdir/dataprocess/knownGene_ucsc_hg19.annot.txt.gz

# Convert format: tab-delimited to bed 
zcat $workdir/dataprocess/knownGene_ucsc_hg19.annot.txt.gz | awk -F$'\t' 'BEGIN{OFS="\t"} {if($2 ~ "chr[0-9]+") print($2,$6,$7,$5)}' | gzip -c > $workdir/dataprocess/knownGene_ucsc_hg19.annot.bed.gz

# 2.2 MiRNA
# Data source: Ensembl BioMart 
wget -O $workdir/rawdata/mirna_targets.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_mirna_target_feature" interface = "default" ><Attribute name = "chromosome_name" /><Attribute name = "chromosome_start" /><Attribute name = "chromosome_end" /><Attribute name = "accession" /></Dataset></Query>'

# Convert format: tab-delimited to bed 
awk -F '\t' 'BEGIN{OFS="\t"} {print "chr"$1,$2-1,$3,"mirna"}' $workdir/rawdata/mirna_targets.txt | gzip -c > $workdir/dataprocess/mirna_targets.bed.gz

# 2.3 Regulatory features
# Data source: Ensembl BioMart 
wget -O $workdir/rawdata/regulatory_features.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_regulatory_feature" interface = "default" ><Attribute name = "chromosome_name" /><Attribute name = "chromosome_start" /><Attribute name = "chromosome_end" /><Attribute name = "regulatory_stable_id" /><Attribute name = "feature_type_name" /></Dataset></Query>'

# Convert format: tab-delimited to bed 
awk -F '\t' 'BEGIN{OFS="\t"} {print "chr"$1,$2-1,$3,"tfbs"}' $workdir/rawdata/regulatory_features.txt | gzip -c > $workdir/dataprocess/regulatory_features.bed.gz

#More details in http://grch37.ensembl.org/biomart/martview/891d411f1bfeb1cd8d2a0f46fe62cff2?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_mirna_target_feature.default.mirna_target_feature.chromosome_name|hsapiens_mirna_target_feature.default.mirna_target_feature.chromosome_start|hsapiens_mirna_target_feature.default.mirna_target_feature.chromosome_end|hsapiens_mirna_target_feature.default.mirna_target_feature.accession&FILTERS=&VISIBLEPANEL=attributepanel

# 2.4 Merge annotation and sort (genes, miRNA and TFBS)
zcat $workdir/dataprocess/knownGene_ucsc_hg19.annot.bed.gz $workdir/dataprocess/mirna_targets.bed.gz $workdir/dataprocess/regulatory_features.bed.gz | sort -k1,1 -k2,2n | gzip -c > $workdir/dataprocess/knownGene_ucsc_hg19.annot.complete.sorted.bed.gz

# 2.5 Sort the genotype file
cat $workdir/dataprocess/chr[0-9]*.bim | awk 'BEGIN{OFS="\t";} {print("chr"$1, $4-1, $4-1+length($5), $2)}' | sort -k1,1 -k2,2n | gzip -c > $workdir/dataprocess/template.1000genomesEAS.sorted.bed.gz

# 2.6 Intersect genotype and annotation
# Intersect 
$workdir/toolkit/bedtools2/bin/bedtools intersect -a $workdir/dataprocess/template.1000genomesEAS.sorted.bed.gz -b $workdir/dataprocess/knownGene_ucsc_hg19.annot.complete.sorted.bed.gz -wa -wb -sorted | gzip -c > $workdir/dataprocess/template.1000genomesEAS.complete_annot_hg19.intersect.txt.gz

# Convert format
python2 $workdir/toolkit/annot2annomat.py $workdir/dataprocess/template.1000genomesEAS.complete_annot_hg19.intersect.txt.gz $workdir/dataprocess/template.1000genomesEAS.sorted.bed.gz $workdir/dataprocess/knownGene_ucsc_hg19.annomat.txt.gz

# 2.7 Prioritize functional annotation 
python2 $workdir/toolkit/uniq_annot.py $workdir/dataprocess/knownGene_ucsc_hg19.annomat.txt.gz $workdir/dataprocess/knownGene_ucsc_hg19.annomat.uniq.txt.gz

# 2.8. Calculate LD r2 including inter-chr LD
mkdir $workdir/dataprocess/schork
for ((i=1;i<23;i++)); do plink --bfile $workdir/dataprocess/chr${i} --r2 inter-chr gz yes-really --ld-window-r2 0.2 --out $workdir/dataprocess/schork/chr${i}.schork.r2; done

# 2.9. Create ld-induced categories.
# This step needs 64G memory. Suggest you run this step on HPC system.
python3 $workdir/toolkit/ld_informed_annot.py $workdir/dataprocess/knownGene_ucsc_hg19.annomat.uniq.txt.gz $workdir/dataprocess/schork/ $workdir/dataprocess/knownGene_ucsc_hg19.annot.ld_informed.txt.gz


# 3.Creat the ref file for pleioFDR (e.g. ref9545380_1kgPhase3eur_LDr2p1.mat)
# Run the following script on Python2.7.
import pandas as pd
import numpy as np
import scipy.io as sio

print('Reading reference...')
ref = pd.read_csv('$workdir/dataprocess/all_chr.bim', delim_whitespace=True, header=None, names='CHR SNP GP BP A1 A2'.split())
ref['A1'] = ref['A1'].str.upper()
ref['A2'] = ref['A2'].str.upper()
snp_to_id = dict([(snp, index) for snp, index in zip(ref['SNP'], ref.index)])

_base_complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
def _complement(seq):
    if any([(b not in _base_complement) for b in seq]): return seq
    return "".join([_base_complement[b] for b in seq])

def _reverse_complement(seq):
    return _complement(seq[::1])

is_ambiguous = [(a1 == _reverse_complement(a2)) for (a1, a2) in zip(ref['A1'],ref['A2'])]

print('Reading allele frequencies...')
frq=pd.concat([pd.read_table('$workdir/dataprocess/chr{}.maf.frq'.format(chri),delim_whitespace=True) for chri in range(1, 23)])

print('Reading LD-weighted annotations...')
intergenic=pd.read_table('$workdir/dataprocess/knownGene_ucsc_hg19.annot.ld_informed.txt.gz', delim_whitespace=True)

ld_window_r2=0.1

ref.to_csv('$workdir/dataprocess/8452254.ref',index=False, sep='\t')

for chri in range(1, 23):
	print('Processing chr{}.r2.ld.gz...'.format(chri))
	df=pd.read_csv('$workdir/dataprocess/chr{}.r2.ld.gz'.format(chri), delim_whitespace=True, usecols=['SNP_A', 'SNP_B', 'R2'])
	df = df[df['R2'] >= ld_window_r2][['SNP_A', 'SNP_B']].copy() # filter based on r2 threshold
	index_A = df['SNP_A'].map(snp_to_id)
	index_B = df['SNP_B'].map(snp_to_id)
	sio.savemat('$workdir/dataprocess/chr{}.r2.ld.mat'.format(chri), {'id1':index_A.values+1, 'id2':index_B.values+1,'chrnumvec':ref['CHR'].values,'posvec':ref['BP'].values,'mafvec':frq['MAF'].values,'is_ambiguous':is_ambiguous,'is_intergenic':intergenic['Intergenic'].values},format='5', do_compression=False, oned_as='column')


# 4.Concatenate data from step3 & convert it to the LD matrix.
# This step needs near 128G memory. Suggest you run this step on HPC system. 
# Run the following script on matlab.
'''
data = load('$workdir/dataprocess/chr1.r2.ld.mat'); id1 = data.id1; id2 = data.id2;
for i=2:22, data = load(sprintf('/$workdir/dataprocess/chr%i.r2.ld.mat', i)) ; id1 = [id1; data.id1]; id2 = [id2; data.id2]; end
nsnp = length(data.chrnumvec)
LDmat = sparse(double(id1),double(id2),true,double(nsnp),double(nsnp));
LDmat = LDmat | speye(double(nsnp));
LDmat = LDmat | (LDmat - LDmat');
data = rmfield(data,{'id1', 'id2'})
data.LDmat = LDmat

save('$workdir/ref8452254_1kgPhase3eas_LDr2p1.mat', '-struct', 'data', '-v7.3')
'''





