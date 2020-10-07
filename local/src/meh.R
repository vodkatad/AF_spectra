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

ma2 <- '/scratch/trcanmed/AF_spectra/dataset/MutationalPattern_bulk/list_fvcf'
ma2_f_comma <- readLines(ma2)
ma2_dir <- '/scratch/trcanmed/AF_spectra/dataset/'
ma2_files <- unlist(strsplit(ma2_f_comma, ','))


pos <- regexpr("CRC\\d+", ma2_files, perl=TRUE)
sample_names_2 <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, ma2_files, pos),'_bulk')

sample_names <- sample_names_2
files <- ma2_files
# starting with mut pat
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

nrun_estimate <- 50
nrun <- 150
seed <- 123456
# we do not add 0.0001 to mut_mat because extract_signatures already do so (other fun used later TODO FIXME?)
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:5, method="brunet", nrun=nrun_estimate, seed=seed)

nsign <- 2 # from the plot optimal n. is 2 considering cophenetic...
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)


names_sign <- paste0("PDO-MA_", seq(1, nsign))
colnames(nmf_res$signatures) <- names_sign
rownames(nmf_res$contribution) <- names_sign

#compare with cosmic
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

#### clevers
# getting list of vcfs
ma2 <- '/scratch/trcanmed/AF_spectra/dataset/MutationalPattern_bulk/list_fvcf'
ma2_f_comma <- readLines(ma2)
ma2_dir <- '/scratch/trcanmed/AF_spectra/dataset/'
ma2_files <- unlist(strsplit(ma2_f_comma, ','))
files_3 <- c("/home/egrassi/P1.merged.vcf.gz","/home/egrassi/P2.merged.vcf.gz","/home/egrassi/P3.merged.vcf.gz")
samples_names_3 <- c('cl1',"cl2","cl3")
pos <- regexpr("CRC\\d+", ma2_files, perl=TRUE)
sample_names_2 <- paste0(mapply(function(x, p) {substr(x,p,p+6)}, ma2_files, pos),'_bulk')

samples_names <- c(sample_names_2, samples_names_3)
files <- c(ma2_files, files_3)
vcfs <- read_vcfs_as_granges(files, samples_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
plot_spectrum(type_occurrences, CT = TRUE, by=samples_names)


nrun_estimate <- 50
nrun <- 150
seed <- 123456
# we do not add 0.0001 to mut_mat because extract_signatures already do so (other fun used later TODO FIXME?)
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:5, method="brunet", nrun=nrun_estimate, seed=seed)

nsign <- 2 # from the plot optimal n. is 2 considering cophenetic...
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)


names_sign <- paste0("PDO-MA_", seq(1, nsign))
colnames(nmf_res$signatures) <- names_sign
rownames(nmf_res$contribution) <- names_sign

#compare with cosmic
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

plot_contribution_heatmap(nmf_res$contribution)
