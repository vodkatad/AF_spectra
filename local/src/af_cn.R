#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
mutfile <- args[1]
maxcn <- as.numeric(args[2])
outputplot <- args[3]
#save.image("p.RData")

muts <- read.table(gzfile(mutfile), header=FALSE, stringsAsFactors=FALSE)
colnames(muts) <- c("chr","b","e","id")
muts$af <- as.numeric(sapply(muts$id, function(x) { strsplit(x, ":")[[1]][7]}))
muts$cn <- sapply(muts$id, function(x) { strsplit(x, ":")[[1]][8]})
muts$cn <- as.numeric(muts$cn)
muts <- muts[muts$cn <= maxcn,]
muts$cn <- factor(muts$cn, levels=seq(1, maxcn))
save.image('pippo.Rdata')
ggplot(data=muts, aes(x=af, y=..count.., color=cn)) + geom_density(position="identity")+theme_bw()
ggsave(outputplot)

