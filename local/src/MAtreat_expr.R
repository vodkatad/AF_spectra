

fp <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/fpkm.tsv.gz'), sep="\t", header=T)

dnds <- read.table('/scratch/trcanmed/AF_spectra/dataset_MAtreats/T_CRC1430/platypus_nobin_00/dnsd_annot.tsv', sep="\t", header=T)

w <- paste0('H_', unique(dnds$gene))
length(w)

fpw <- fp[w, grepl('CRC1430', colnames(fp))]

summary(t(fpw))

# 10 su 34 non espressi.