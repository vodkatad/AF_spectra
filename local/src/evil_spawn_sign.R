library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(pheatmap)
vcfs <- c('/scratch/trcanmed/AF_spectra/dataset/CRC1502_clones_all/platypus_nobin/private_9/0001.vcf','/scratch/trcanmed/AF_spectra/dataset/CRC1502_clones_all/platypus_nobin/private_9/0000.vcf')
sample_names <- c('private_evil9','other_2nd_clones')
vcf_files <- vcfs
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
cosmic <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "") # ???
sp_url <- paste(cosmic, sep = "")
ref_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), ref_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe>
ref_signatures = ref_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names>
row.names(ref_signatures) = ref_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
ref_signatures = as.matrix(ref_signatures[,4:33])
##
## signal TODO  new version of package has it already inside?
## other cosmic version
###
ff <- fit_to_signatures(mut_mat, ref_signatures)
ff$contribution
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
xlab("") +theme_bw()+ theme(axis.text.y = element_text(size=10), axis.title.x=element_text(size=15))+
geom_hline(aes(yintercept=.95)) +scale_fill_brewer(palette = "Dark2")
data <- as.matrix(ff$contribution)
#select <- which(rowSums(data) > 10)
data = t(data)
data = data/rowSums(data)
#select <- which(colSums(data)>1)
#data <- data[,select]
annot_rows <- data.frame(row.names=rownames(data))
head(annot_rows)
head(data)
rownames(data)
col <- brewer.pal(3,'Dark2')
pheatmap(data, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE,  color=brewer.pal(9, 'PuBu'), cluster_rows=FALSE)
select <- which(colSums(data)!=0)
data <- data[,select]
pheatmap(data, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE,  color=brewer.pal(9, 'PuBu'), cluster_rows=FALSE)

