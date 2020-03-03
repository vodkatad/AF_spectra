#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
mutfile <- args[1]
maxcn <- as.numeric(args[2])
outputplot <- args[3]
muts <- read.table(gzfile(mutfile), header=FALSE, stringsAsFactors=FALSE)
colnames(muts) <- c("chr","b","e","id")
muts$af <- as.numeric(sapply(muts$id, function(x) { strsplit(x, ":")[[1]][7]}))
muts$cn <- sapply(muts$id, function(x) { strsplit(x, ":")[[1]][8]})
muts$cn <- as.numeric(muts$cn)
muts <- muts[muts$cn <= maxcn,]
muts$cn <- factor(muts$cn, levels=seq(1, maxcn))
muts$naf <- muts$af * muts$cn
cat(dim(muts))
cat("\n")
muts <- muts[muts$naf <=1,]
cat(dim(muts))
cat("\n")
ggplot(data=muts, aes(x=naf)) + geom_density()+theme_bw()
ggsave(outputplot)

save.image('af_normalize.Rdata')