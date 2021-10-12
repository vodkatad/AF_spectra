#!/usr/bin/env Rscript
library(UpSetR)
args <- commandArgs(trailingOnly = T)
input <- args[1]
output <- args[2]
save.image('pippo.Rdata')
d <- read.table(input, sep="\t", header=FALSE)
colnames(d) <- c('mut','sample')
listmut <- lapply(unique(d$sample), function(x) {d[d$sample==x,'mut']})
names(listmut)  <- unique(d$sample)
svg(output)
upset(fromList(listmut), order.by = "freq", sets=sort(names(listmut)), keep.order=TRUE, text.scale=0.5)
graphics.off()