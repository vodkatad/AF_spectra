#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]

data <- read.table(gzfile(infile), sep="\t", header=TRUE, quote="")
rownames(data) <- paste0(data[,'CHROM'],":", data[,'POS'], ":", data[,'REF'], ":", data[,'ALT'])
data[,'CHROM'] <- NULL
data[,'POS'] <- NULL
data[,'REF'] <- NULL
data[,'ALT'] <- NULL

not_private <- rowSums(data) > 1

data <- data[not_private,]
data$fake <- rep(0, nrow(data))
write.table(data.frame("ID"=rownames(data),data), file = gzfile(outfile), sep="\t", quote=FALSE, row.names=FALSE)