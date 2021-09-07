print('This script would take 6 arguments, which are:')
print('Conditional clump 0.01 lead csv, respective FUMA snps.tx, std sumstats for TRAIT1')
print('std sumstats for TRAIT2, TRAIT1 name, TRAIT2 name')
print('to produce cond_0.01_TRAIT1_vs_TRAIT2.csv')
args <- commandArgs(TRUE)
#lead1=read.csv('fdr_out5/SCZ_BMI_condfdr/cond.result.clump_0.01.lead.csv',sep = '\t')
lead1=read.csv(args[1],sep = '\t')
lead2=subset(lead1, is_locus_lead == "True")
lead3=lead2[order(lead2$locusnum, lead2$FDR),]
lead4 = lead3[!duplicated(lead3$locusnum),]
lead5=subset(lead4, select = -c(is_locus_lead) )
#snp1=read.csv('fuma_result_fin1/FUMA_cond_scz_bmi_lead/snps.txt',sep = '\t')
snp1=read.csv(args[2],sep = '\t')
leadsnp1=merge(lead5,snp1,by.x = 'LEAD_SNP',by.y = 'rsID')
leadsnp2=leadsnp1[,c(2,3,1,4,5,6,7,12,11,18,19,20,21,22,23,24)]
#sm1=read.table('sumstat-std/PGC_SCZ_0518_EUR.sumstats.gz',sep = '\t',header = T)
#sm2=read.table('sumstat-std/GIANT_BMI_2018_UKB_v2.sumstats.gz',sep = '\t',header = T)
sm1=read.table(args[3],sep = '\t',header = T)
sm2=read.table(args[4],sep = '\t',header = T)
trait1=args[5]
trait2=args[6]
#trait1='SCZ'
#trait2='BMI'
#st_cols=list('PVAL','Z','OR','BETA','SE')

if ('PVAL' %in% colnames(sm1) & ('PVAL' %in% colnames(sm2)))
  leadsnp2[paste0('PVAL_in_', trait1)]=sm1$PVAL[match(leadsnp2$LEAD_SNP,sm1$SNP)]
if ('Z' %in% colnames(sm1) & ('Z' %in% colnames(sm2)))
  leadsnp2[paste0('Z_in_', trait1)]=sm1$Z[match(leadsnp2$LEAD_SNP,sm1$SNP)]
if ('OR' %in% colnames(sm1) & ('OR' %in% colnames(sm2)))
  leadsnp2[paste0('OR_in_', trait1)]=sm1$OR[match(leadsnp2$LEAD_SNP,sm1$SNP)]
if ('BETA' %in% colnames(sm1) & ('BETA' %in% colnames(sm2)))
  leadsnp2[paste0('BETA_in_', trait1)]=sm1$BETA[match(leadsnp2$LEAD_SNP,sm1$SNP)]
if ('SE' %in% colnames(sm1) & ('SE' %in% colnames(sm2)))
  leadsnp2[paste0('SE_in_', trait1)]=sm1$SE[match(leadsnp2$LEAD_SNP,sm1$SNP)]

if ('PVAL' %in% colnames(sm1) & ('PVAL' %in% colnames(sm2)))
  leadsnp2[paste0('PVAL_in_', trait2)]=sm2$PVAL[match(leadsnp2$LEAD_SNP,sm2$SNP)]
if ('Z' %in% colnames(sm1) & ('Z' %in% colnames(sm2)))
  leadsnp2[paste0('Z_in_', trait2)]=sm2$Z[match(leadsnp2$LEAD_SNP,sm2$SNP)]
if ('OR' %in% colnames(sm1) & ('OR' %in% colnames(sm2)))
  leadsnp2[paste0('OR_in_', trait2)]=sm2$OR[match(leadsnp2$LEAD_SNP,sm2$SNP)]
if ('BETA' %in% colnames(sm1) & ('BETA' %in% colnames(sm2)))
  leadsnp2[paste0('BETA_in_', trait2)]=sm2$BETA[match(leadsnp2$LEAD_SNP,sm2$SNP)]
if ('SE' %in% colnames(sm1) & ('SE' %in% colnames(sm2)))
  leadsnp2[paste0('SE_in_', trait2)]=sm2$SE[match(leadsnp2$LEAD_SNP,sm2$SNP)]

names(leadsnp2)[8]='A1'
names(leadsnp2)[9]='A2'

leadsnp3=subset(leadsnp2, FDR < 0.01)
leadsnp4=leadsnp3[order(leadsnp3$locusnum),]
write.csv(leadsnp4, file = paste0('condFDR_0.01_',trait1,'_vs_',trait2, ".csv"), sep='\t', row.names = FALSE)
