#!/usr/bin/env Rscript

inputs <- snakemake@input

load_polish <- function(filename) {
    first <- read.table(filename, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    colnames(first) <- c('chr', 'b', 'e', 'id')
    first$idnoAF <- sapply(first$id, function(x) {y<-strsplit(x, ':')[[1]]; return(paste0(y[c(1,2,3,4)], collapse=":"))})
    first
}

save.image('wtf.Rdata')
sample_mut <- data.frame(sample=c(), mut=c())
for (i in seq(1, length(inputs))) {
    this_d <- load_polish(inputs[[i]])
    sample_mut <- rbind(sample_mut, data.frame(sample=rep(i, nrow(this_d)), mut=this_d$idnoAF))
}

tabledf <- as.data.frame(table(sample_mut$mut))
colnames(tabledf) <- c('mut', 'n')

res <- tabledf[tabledf$n == length(inputs),, drop=FALSE] 
res$nn <- rep('common', length(res))
res$n <- NULL

# efficiency and good code at their peaks
res2 <- tabledf[tabledf$n == 1,, drop=FALSE]
m <- merge(res2, sample_mut, by="mut")
m$nn <- paste0('c', m$sample)
m$n <- NULL
m$sample <- NULL
res <- rbind(res, m)
write.table(res, gzfile(snakemake@output[['res']]), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)