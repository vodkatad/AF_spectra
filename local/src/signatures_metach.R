library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)
library(pheatmap)

# getting list of vcf
ma_files <- c('n1.vcf.gz', 'n2.vcf.gz')
ma_dir <- "/scratch/trcanmed/AF_spectra/dataset/CRC1078/tree/"
files <- c(paste0(ma_dir, ma_files))

sample_names <- c("T0_private_CRC1078","T0_private_CRC1307")

# starting with mut pat
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

### vs signal
signal <- '/home/egrassi/signal_colorectal.txt'
crc_signatures = read.table(signal, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), crc_signatures$Substitution.Type)
# Reorder cancer signatures dataframe>
crc_signatures = crc_signatures[as.vector(new_order),]
crc_signatures = as.matrix(crc_signatures[,2:(ncol(crc_signatures)-1)])
# Add trinucletiode changes names as row.names>
row.names(crc_signatures) = crc_signatures$Substitution.Type
ff <- fit_to_signatures(mut_mat, crc_signatures)
ff$contribution
which(rowSums(ff$contribution) > 10)
select <- which(rowSums(ff$contribution) > 10)
plot_contribution_heatmap(ff$contribution,cluster_samples = TRUE,method = "complete")



cos_sim_ori_rec <- cos_sim_matrix(mut_mat, ff$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$s = row.names(cos_sim_ori_rec)
cos_sim_ori_rec$sample <- unlist(lapply(strsplit(cos_sim_ori_rec$s,"_"), function(x) x[length(x)]))
cos_sim_ori_rec$n <- seq(1, nrow(cos_sim_ori_rec))
ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=reorder(s, n), fill=sample)) +
  geom_bar(stat="identity") +
  coord_flip(ylim=c(0.7,1)) +
  ylab("Cosine similarity\n original VS reconstructed") + 
  xlab("") +theme_bw()+ theme(axis.text.y = element_text(size=15), axis.title.x=element_text(size=15))+
  geom_hline(aes(yintercept=.95)) +scale_fill_brewer(palette = "Dark2")