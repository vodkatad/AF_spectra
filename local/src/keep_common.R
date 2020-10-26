#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
time_to_keep <- args[3]

data <- read.table(gzfile(infile), sep="\t", header=TRUE, quote="")
rownames(data) <- paste0(data[,'CHROM'],":", data[,'POS'], ":", data[,'REF'], ":", data[,'ALT'])
data[,'CHROM'] <- NULL
data[,'POS'] <- NULL
data[,'REF'] <- NULL
data[,'ALT'] <- NULL
save.image('pipo.Rdata')
# - became . when loading with read.table
time <- sapply(strsplit(colnames(data), '.', fixed=TRUE), function(x) { x[3] })
data <- data[, time==time_to_keep]
if (nrow(data) != 0) {
    common <- rowSums(data) == ncol(data)
    data <- data[common,]
    write.table(data.frame("ID"=rownames(data),data), file = gzfile(outfile), sep="\t", quote=FALSE, row.names=FALSE)
} else {
    write.table(data.frame("wrongtimeselected"=c("NA")), file = gzfile(outfile), sep="\t", quote=FALSE, row.names=FALSE)
}