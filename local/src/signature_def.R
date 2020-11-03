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


### cosmic single clones TODO + bulk
# TODO rules then with cosmic 6
# TODo SigProfilerExtractor

#egrassi@godot:/scratch/trcanmed/AF_spectra/dataset/AF_spectra$ tr "," "\n" < gained_clones | grep M | tr "\n" "," > gained_clones_vivo
#egrassi@godot:/scratch/trcanmed/AF_spectra/dataset/AF_spectra$ tr "," "\n" < gained_clones | grep -v M | tr "\n" "," > gained_clones_vitro

### sigle clones vitro analyses
f <- '/scratch/trcanmed/AF_spectra/dataset/AF_spectra/gained_clones_vitro'
f_comma <- readLines(f)
files <- unlist(strsplit(f_comma, ','))
# assembling samples names
pos <- regexpr("CRC\\d+\\-", files, perl=TRUE)
sample_names <- paste0(mapply(function(x, p) {substr(x,p,p+13)}, files, pos), "_MA_vitro")
###  bulk
bulk <- '/scratch/trcanmed/AF_spectra/dataset/MutationalPattern_bulk/list_vcf'
bulk_f_comma <- readLines(bulk)
bulk_dir <- '/scratch/trcanmed/AF_spectra/dataset/'
bulk_files <- unlist(strsplit(bulk_f_comma, ','))
pos <- regexpr("CRC\\d+", bulk_files, perl=TRUE)
sample_names_bulk <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, bulk_files, pos),'_bulk')

files <- c(files, bulk_files)
sample_names <- c(sample_names, sample_names_bulk)

fv <- '/scratch/trcanmed/AF_spectra/dataset/AF_spectra/gained_clones_vivo'
fv_comma <- readLines(fv)
files_vivo <- unlist(strsplit(fv_comma, ','))
# assembling samples names
pos <- regexpr("CRC\\d+\\-", files_vivo, perl=TRUE)
sample_names_vivo <- paste0(mapply(function(x, p) {substr(x,p,p+14)}, files_vivo, pos), "_MA_vivo")

files <- c(files, files_vivo)
sample_names <- c(sample_names, sample_names_vivo)

# starting with mut pat
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

# cosmic signatures
cosmic <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "") # ???
sp_url <- paste(cosmic, sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe> 
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names>
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])


ff <- fit_to_signatures(mut_mat, cancer_signatures)
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
  xlab("") +theme_bw()+ theme(axis.text.y = element_text(size=10), axis.title.x=element_text(size=15))+
  geom_hline(aes(yintercept=.95)) +scale_fill_brewer(palette = "Dark2")


data <- as.matrix(ff$contribution)
#select <- which(rowSums(data) > 10)
data = t(data)
data = data/rowSums(data)
select <- which(colSums(data)>1)
data <- data[,select]
annot_rows <- data.frame(row.names=rownames(data))
#annot_rows[grepl('M',rownames(annot_rows)), 'sample'] <- 'vivo'
annot_rows$sample <-  as.factor(unlist(lapply(strsplit(rownames(annot_rows),'_'), function(x){ x[length(x)] })))
annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] }))
annot_rows$model <- as.factor(unlist(lapply(strsplit(annot_rows$model,'_'), function(x){ x[1] })))
#annot_cols <- data.frame(signature=c('Clock','Clock','ROS','smoke','MMR','MMR','NA','NA'), row.names=colnames(data))
col <- brewer.pal(3,'Dark2')
names(col) <- levels(annot_rows$sample)

cbPalette2 <- c("#ff5733", 
                #"#9d01fc"
                "#f607b9",
                "#155d00",
                "#77a003",
                "#0829fc"
)
names(cbPalette2) <- levels(annot_rows$model)
annot_colors <- list(sample=col, model=cbPalette2)
pheatmap(data, fontsize_row = 9, fontsize_col=13, show_colnames = TRUE,  cluster_cols=FALSE, angle_col = 315,  annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(50, 'PuBu'))

