#!/usr/bin/env Rscript

library(dndscv)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
image <- args[3]
rda <- args[4]

getdnds <- function(file) {
	data <- read.table(file, sep="\t", header=FALSE)
	colnames(data) <- c('sampleID','chr', 'pos', 'ref', 'alt')
	return(dndscv(data, refdb=rda, cv=NULL))
	#return(dndscv(data, refdb='hg19', cv='hg19'))
}

res <- getdnds(infile)
write.table(res$globaldnds, file = outfile, sep="\t", quote=FALSE)
save.image(image)