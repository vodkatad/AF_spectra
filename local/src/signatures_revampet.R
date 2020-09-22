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
tcga <- '/scratch/trcanmed/pdxopedia/dataset/tcga/list_vcf'
ma <- '/scratch/trcanmed/AF_spectra/dataset/list_vcf_vitro'
ma_f_comma <- readLines(ma)
tcga_f_comma <- readLines(tcga)
tcga_dir <- '/scratch/trcanmed/pdxopedia/dataset/tcga/'
ma_dir <- '/scratch/trcanmed/AF_spectra/dataset/'
ma_files <- unlist(strsplit(ma_f_comma, ','))
tcga_files <- unlist(strsplit(tcga_f_comma, ','))
files <- c(paste0(ma_dir, ma_files), paste0(tcga_dir, tcga_files))

# assembling samples names
pos <- regexpr("CRC\\d+", ma_files, perl=TRUE)
sample_names <- mapply(function(x, p) {substr(x,p,p+6)}, ma_files, pos)
sample_names <- c(sample_names, paste0('tcga', seq(1, length(tcga_files))))

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

