#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
colors <- args[3]

d <- read.table(infile, sep="\t", quote="", header=FALSE)
colnames(d) <- c('model','estimate', 'upper', 'lower')
d <- d[d$model != "overall",]

cbPalette2 <- unlist(strsplit(colors, ','))

ggplot(d, aes(x=model, y=estimate, color=model)) +  geom_point(stat="identity", size=3) +
    geom_errorbar(aes(ymin=lower, ymax=upper, x=model, width=0.1, color=model), size=0.6)+theme_bw()+ggtitle('dN/dS nonsyn+truncating')+ylab('ML estimate')+
  theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), legend.position="none",axis.title.y=element_text(size=15))+scale_color_manual(values=cbPalette2)

ggsave(outfile)