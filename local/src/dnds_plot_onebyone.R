#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = T)
infile <- args[1] # need to merge with ','
outfile <- args[2]

TODO

models <- c('CRC0282','CRC0327','CRC1078', 'CRC1307','CRC1502')

plotdndssingle <- function(model) {
  d <- read.table(paste0('~/work/evol/MA/',model,'_dnds.tsv'), sep="\t", header=T)

  ggplot(d, aes(x=name, y=mle, color=name)) +  geom_point(stat="identity", size=3) +
  geom_errorbar(aes(ymin=cilow, ymax=cihigh, x=name, width=0.1, color=name), size=0.6)+theme_bw()+ylab('ML estimate')+
  theme(axis.text.x = element_text(size=15), legend.position="none")+ylim(-1,4.8)+ggtitle('')
}

pl <- lapply(models, plotdndssingle)

library(ggpubr)
ggarrange(pl[[1]],pl[[2]],pl[[3]],pl[[4]],pl[[5]],
                    labels = models,
                   ncol = 2, nrow=3)
