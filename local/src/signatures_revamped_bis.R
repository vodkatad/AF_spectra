# us vs in vivo
library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)
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
sample_names <- mapply(function(x, p) {substr(x,p,p+6)}, ma_files, pos)
pos <- regexpr("CRC\\d+", ma2_files, perl=TRUE)
sample_names_2 <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, ma2_files, pos),'_bulk')
sample_names <- c(sample_names, sample_names_2)
# starting with mut pat
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
# NMF de novo signature extraction
nrun_estimate <- 50
nrun <- 150
seed <- 123456
# we do not add 0.0001 to mut_mat because extract_signatures already do so (other fun used later TODO FIXME?)
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:5, method="brunet", nrun=nrun_estimate, seed=seed)
nsign <- 4 # from the plot optimal n. is 2 considering cophenetic...
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)
names_sign <- paste0("PDO-MA_", seq(1, nsign))
colnames(nmf_res$signatures) <- names_sign
rownames(nmf_res$contribution) <- names_sign
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
us <- nmf_res$signatures
colnames(us) <- names_sign
hclust_cosmic_us = cluster_signatures(cbind(us, cancer_signatures), method = "average")
plot(hclust_cosmic_us)
ci <- c('MA1_cosmic8','MA2_cosmic20','MA3_cosmic9','MA3_cosmic1-6')
colnames(nmf_res$signatures) <- ci
rownames(nmf_res$contribution) <- ci
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
sample_names <- mapply(function(x, p) {substr(x,p,p+6)}, ma_files, pos)
pos <- regexpr("CRC\\d+", ma2_files, perl=TRUE)
sample_names_2 <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, ma2_files, pos),'_bulk')
sample_names <- c(sample_names, sample_names_2)
fs <- c('/scratch/trcanmed/AF_spectra/dataset/CRC0282/platypus_nobin/vivo.merged.vcf.gz','/scratch/trcanmed/AF_spectra/dataset/CRC1307_platypus_nobin/vivo.merged.vcf.gz')
vivo_files <- unlist(strsplit(fs, ','))
pos <- regexpr("CRC\\d+", vivo_files, perl=TRUE)
sample_names_3 <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, vivo_files, pos),'_vivo')
files <- c(files, vivo_files)
sample_names <- c(sample_names, sample_names_3)
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
ff <- fit_to_signatures(mut_mat,  nmf_res$signatures)
ff$contribution
plot_contribution_heatmap(ff$contribution,cluster_samples = TRUE,method = "complete")

# re train the signatures using also in vivo data
# same phenomenon of less 8 in CRC1307_vivo but not huge

############ lost
# getting list of vcfs
ma <- '/scratch/trcanmed/AF_spectra/dataset/list_vcf_vitroloss'
ma_f_comma <- readLines(ma)
ma_dir <- '/scratch/trcanmed/AF_spectra/dataset/'
ma_files <- unlist(strsplit(ma_f_comma, ','))
files <- c(paste0(ma_dir, ma_files))
# assembling samples names
pos <- regexpr("CRC\\d+", ma_files, perl=TRUE)
sample_names <- mapply(function(x, p) {substr(x,p,p+6)}, ma_files, pos)
# starting with mut pat
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat_lost <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)


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
sample_names <- mapply(function(x, p) {substr(x,p,p+6)}, ma_files, pos)
pos <- regexpr("CRC\\d+", ma2_files, perl=TRUE)
sample_names_2 <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, ma2_files, pos),'_bulk')
sample_names <- c(sample_names, sample_names_2)
# starting with mut pat
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)


signal <- '/home/egrassi/signal_colorectal.txt'
crc_signatures = read.table(signal, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), crc_signatures$Substitution.Type)
# Reorder cancer signatures dataframe> 
crc_signatures = crc_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names>
row.names(crc_signatures) = crc_signatures$Substitution.Type
# Keep only 96 contributions of the signatures in matrix
crc_signatures = as.matrix(crc_signatures[,2:(ncol(crc_signatures)-1)])
# manually checked that the last one was never there, POLE.

ff <- fit_to_signatures(mut_mat, crc_signatures)

  garb <- sapply(sample_names, function(x) { 
    #pdf(paste0("cosmic_",x,".signature_evaluation.pdf"))
    pos <- which(sample_names==x)
    ggsave(file=paste0("cosmic_",x,".signature_evaluation.png"), plot=plot_compare_profiles(mut_mat[,pos], ff$reconstructed[,pos],profile_names = c("Original", "Reconstructed"), condensed = TRUE))
    #graphics.off()
  })
