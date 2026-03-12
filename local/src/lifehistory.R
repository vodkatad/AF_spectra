library(ggplot2)
s4c <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/S4C.tsv', sep="\t", header=T)
s4d <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/S4D.tsv', sep="\t", header=T)


table(s4c$tumour_type)
s4c_coread <- s4c[s4c$tumour_type == 'ColoRect-AdenoCa',]
s4d_coread <- s4d[s4d$tumour_type == 'ColoRect-AdenoCa',]


ggplot(data=s4c_coread, aes(x=signature, y=log2fc_earlyLate))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(height=0, size=1)+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept=0, color='red')


ggplot(data=s4d_coread, aes(x=signature, y=log2fc_clonalSubclonal))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(height=0, size=1)+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept=0, color='red')

ggplot(data=s4d, aes(x=signature, y=log2fc_clonalSubclonal))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(height=0, size=1)+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept=0, color='red')


ggplot(data=s4c, aes(x=signature, y=log2fc_earlyLate))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(height=0, size=1)+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept=0, color='red')



sign <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/PCAWG_sub_signatures_in_samples_beta2.20170320.donor', sep='\t', header=T)
tt <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/pcawg_specimen_histology_August2016_v9_donor', sep='\t', header=T)

coad_pz <- tt[tt$histology_abbreviation== 'ColoRect-AdenoCA',]

sign_coad <- sign[sign$donor %in% coad_pz$icgc_specimen_id,]

summary(sign_coad$Signature_8)
#> table(sign$donor %in% substr(coad_pz$icgc_specimen_id, 0, 5))

#FALSE  TRUE 
#2707     1 

#???

# https://pcawg-hub.s3.us-east-1.amazonaws.com/download/October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS
# hg19
