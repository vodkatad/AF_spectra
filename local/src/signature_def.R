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

# getting list of vcfs
ma <- '/scratch/trcanmed/AF_spectra/dataset/list_vcf_vitro'
ma_f_comma <- readLines(ma)
ma_dir <- '/scratch/trcanmed/AF_spectra/dataset/'
ma_files <- unlist(strsplit(ma_f_comma, ','))
files <- c(paste0(ma_dir, ma_files))

ma2 <- '/scratch/trcanmed/AF_spectra/dataset/MutationalPattern_bulk/list_vcf'
ma2_f_comma <- readLines(ma2)
ma2_dir <- '/scratch/trcanmed/AF_spectra/dataset/'
ma2_files <- unlist(strsplit(ma2_f_comma, ','))

files <- c(paste0(ma_dir, ma_files), ma2_files)
# assembling samples names
pos <- regexpr("CRC\\d+", ma_files, perl=TRUE)
sample_names <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, ma_files, pos), '_MA_vitro')

pos <- regexpr("CRC\\d+", ma2_files, perl=TRUE)
sample_names_2 <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, ma2_files, pos),'_bulk')

sample_names <- c(sample_names, sample_names_2)

fs <- c('/scratch/trcanmed/AF_spectra/dataset/CRC0282/platypus_nobin/vivo.merged.vcf.gz','/scratch/trcanmed/AF_spectra/dataset/CRC1307_platypus_nobin/vivo.merged.vcf.gz')
pos <- regexpr("CRC\\d+", fs, perl=TRUE)
sample_names3 <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, fs, pos), '_MA_vivo')


files <- c(files, fs)
sample_names <- c(sample_names, sample_names3)

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
ggsave('bulk_merged_cosinesim.svg')

data <- as.matrix(ff$contribution)
data = t(data)
data = data/rowSums(data)
annot_rows <- data.frame(sample= unlist(lapply(strsplit(rownames(data),"_"), function(x) x[length(x)])), row.names=rownames(data))
#annot_cols <- data.frame(signature=c('Clock','Clock','ROS','smoke','MMR','MMR','NA','NA'), row.names=colnames(data))
col <- brewer.pal(3,'Dark2')
names(col) <- levels(annot_rows$sample)
annot_colors <- list(sample=col)
pheatmap(data, fontsize_row = 13, fontsize_col=13, show_colnames = TRUE,  cluster_cols=FALSE, angle_col = 315,  annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(50, 'PuBu'))
