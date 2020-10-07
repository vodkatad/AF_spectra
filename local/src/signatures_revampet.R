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
#nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun, single_core=TRUE)


my_extract_signatures <- function(mut_matrix, rank, nrun = 200, nmf_type = c("regular", "variational_bayes"), single_core = FALSE) {
  # Match argument
  nmf_type <- match.arg(nmf_type)
  
  # Add a small pseudocount to avoid features with zero counts.
  mut_matrix <- as.matrix(mut_matrix) + 0.0001
  
  # Make sure the rank_range is valid.
  if (!(rank > 0 & rank == round(rank))) {
    stop("Rank should be a positive integer", call. = FALSE)
  }
  
  if (ncol(mut_matrix) < max(rank)) {
    stop(paste0(
      "The rank should be smaller than the number of ",
      "samples in the input matrix."
    ), call. = FALSE)
  }
  
  if (nmf_type == "regular") {
    # Calculate NMF
    if (single_core){
      res <- NMF::nmf(mut_matrix, rank = rank, method = "brunet", nrun = nrun, seed = 123456, .opt = "v-p")
    } else{
      res <- NMF::nmf(mut_matrix, rank = rank, method = "brunet", nrun = nrun, seed = 123456)
    }
    # Find signatures and contribution of signatures
    signatures <- NMF::basis(res)
    contribution <- NMF::coef(res)
  } else {
    if (!requireNamespace("ccfindR", quietly = TRUE)) {
      stop(paste0(
        "Package 'ccfindR' is needed for variational_bayes to work. ",
        "Please either install it or use the regular NMF."
      ), call. = FALSE)
    }
    sc <- ccfindR::scNMFSet(count = mut_matrix)
    res <- ccfindR::vb_factorize(sc, ranks = rank, nrun = nrun, progress.bar = FALSE, verbose = 0)
    # estimate = ccfindR::vb_factorize(sc, ranks = 2:7, nrun = nrun, progress.bar = FALSE, verbose = 0)
    # plot(estimate)
    # optimal_rank(sb)
    signatures <- ccfindR::basis(res)[[1]]
    contribution <- ccfindR::coeff(res)[[1]]
  }
  # Reconstruct mutation matrix
  reconstructed <- signatures %*% contribution
  return(list(
    signatures = signatures,
    contribution = contribution,
    reconstructed = reconstructed
  ))
}

nmf_res <- my_extract_signatures(mut_mat, rank = nsign, nrun = nrun, single_core=TRUE)


## TODO match to these signatures our data.
nn <- nmf_res$contribution
nn <- nn[, grepl('CRC',colnames(nn))]
plot_contribution_heatmap(nn)

# todo on gained