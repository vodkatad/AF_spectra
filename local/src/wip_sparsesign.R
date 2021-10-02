library(SparseSignatures)
library("BSgenome.Hsapiens.UCSC.hg38")
setwd('/scratch/trcanmed/AF_spectra/dataset')

inf <- 'MSS_signin.tsv'

bsg <- BSgenome.Hsapiens.UCSC.hg38
data(mutation_categories)
head(mutation_categories)

data <- read.table(inf, sep="\t", header=T)
data$chrom <- paste0('chr', data$chrom) # cause we are using UCSC
imported_data <- import.counts.data(input=data,bsg=bsg,mutation_categories=mutation_categories)