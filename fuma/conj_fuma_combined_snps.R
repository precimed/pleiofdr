print('This script would take 6 arguments, which are:')
print('Conjuntional clump 0.05 snps csv, respective FUMA snps.tx, std sumstats for TRAIT1')
print('std sumstats for TRAIT2, TRAIT1 name, TRAIT2 name')
print('to produce conj_0.05_TRAIT1_vs_TRAIT2_snps.csv')
library(data.table)
args <- commandArgs(TRUE)
#lead1=read.csv('fuma_input/conj.md_bmi.snps.csv',sep = '\t')
lead1=fread(args[1],sep = '\t',stringsAsFactors=F,data.table=F)
#lead2=subset(lead1, is_locus_lead == "True")
lead3=lead1[order(lead1$locusnum, lead1$FDR),]
lead4 = lead3[!duplicated(lead3$CAND_SNP),]
#lead5=subset(lead4, select = -c(is_locus_lead) )
lead5=lead4
#snp1=read.csv('fuma_output_fin1/FUMA_conj_bmi_snps_job139939/snps.txt',sep = '\t')
snp1=fread(args[2],sep = '\t',stringsAsFactors=F,data.table=F)
leadsnp1=merge(lead5,snp1,by.x = 'CAND_SNP',by.y = 'rsID')
leadsnp2=leadsnp1[,c(2,3,4,5,1,6,7,8,9,10,11,20:31)]
#sm1=read.table('sumstat-std/PGC_SCZ_0518_EUR.sumstats.gz',sep = '\t',header = T)
#sm2=read.table('sumstat-std/GIANT_BMI_2018_UKB_v2.sumstats.gz',sep = '\t',header = T)
sm1=fread(args[3],sep = '\t',header = T,stringsAsFactors=F,data.table=F)
sm2=fread(args[4],sep = '\t',header = T,stringsAsFactors=F,data.table=F)
trait1=args[5]
trait2=args[6]
#trait1='SCZ'
#trait2='BMI'
#st_cols=list('PVAL','Z','OR','BETA','SE')

if ('PVAL' %in% colnames(sm1) & ('PVAL' %in% colnames(sm2)))
  leadsnp2[paste0('PVAL_in_', trait1)]=sm1$PVAL[match(leadsnp2$CAND_SNP,sm1$SNP)]
if ('Z' %in% colnames(sm1) & ('Z' %in% colnames(sm2)))
  leadsnp2[paste0('Z_in_', trait1)]=sm1$Z[match(leadsnp2$CAND_SNP,sm1$SNP)]
if ('OR' %in% colnames(sm1) & ('OR' %in% colnames(sm2)))
  leadsnp2[paste0('OR_in_', trait1)]=sm1$OR[match(leadsnp2$CAND_SNP,sm1$SNP)]
if ('BETA' %in% colnames(sm1) & ('BETA' %in% colnames(sm2)))
  leadsnp2[paste0('BETA_in_', trait1)]=sm1$BETA[match(leadsnp2$CAND_SNP,sm1$SNP)]
if ('SE' %in% colnames(sm1) & ('SE' %in% colnames(sm2)))
  leadsnp2[paste0('SE_in_', trait1)]=sm1$SE[match(leadsnp2$CAND_SNP,sm1$SNP)]

if ('PVAL' %in% colnames(sm1) & ('PVAL' %in% colnames(sm2)))
  leadsnp2[paste0('PVAL_in_', trait2)]=sm2$PVAL[match(leadsnp2$CAND_SNP,sm2$SNP)]
if ('Z' %in% colnames(sm1) & ('Z' %in% colnames(sm2)))
  leadsnp2[paste0('Z_in_', trait2)]=sm2$Z[match(leadsnp2$CAND_SNP,sm2$SNP)]
if ('OR' %in% colnames(sm1) & ('OR' %in% colnames(sm2)))
  leadsnp2[paste0('OR_in_', trait2)]=sm2$OR[match(leadsnp2$CAND_SNP,sm2$SNP)]
if ('BETA' %in% colnames(sm1) & ('BETA' %in% colnames(sm2)))
  leadsnp2[paste0('BETA_in_', trait2)]=sm2$BETA[match(leadsnp2$CAND_SNP,sm2$SNP)]
if ('SE' %in% colnames(sm1) & ('SE' %in% colnames(sm2)))
  leadsnp2[paste0('SE_in_', trait2)]=sm2$SE[match(leadsnp2$CAND_SNP,sm2$SNP)]

leadsnp3=subset(leadsnp2, FDR < 0.05 & R2 > 0.6)
leadsnp4=leadsnp3[order(leadsnp3$locusnum),]


leadsnp4[paste0('A1_in_',trait1)]=sm1$A1[match(leadsnp4$CAND_SNP,sm1$SNP)]
leadsnp4[paste0('A1_in_',trait2)]=sm2$A1[match(leadsnp4$CAND_SNP,sm2$SNP)]


leadsnp4[paste0('Z_recalculated_in_',trait1)]=ifelse(leadsnp4$A1 == leadsnp4[[paste0('A1_in_',trait1)]], 
                                                     leadsnp4[[paste0('Z_in_',trait1)]], (-1)*leadsnp4[[paste0('Z_in_',trait1)]])
leadsnp4[paste0('Z_recalculated_in_',trait2)]=ifelse(leadsnp4$A1 == leadsnp4[[paste0('A1_in_',trait2)]], 
                                                     leadsnp4[[paste0('Z_in_',trait2)]], (-1)*leadsnp4[[paste0('Z_in_',trait2)]])


leadsnp4$ConcordEffect=ifelse(leadsnp4[[paste0('Z_recalculated_in_',trait1)]] > 0 & leadsnp4[[paste0('Z_recalculated_in_',trait2)]] > 0,'Yes',
                              ifelse(leadsnp4[[paste0('Z_recalculated_in_',trait1)]] < 0 & leadsnp4[[paste0('Z_recalculated_in_',trait2)]] < 0,
                                     'Yes','No'))
write.table(leadsnp4, file = paste0('conjFDR_0.05_',trait1,'_vs_',trait2,'_snps', ".csv"), sep='\t', row.names = F,quote=F)
