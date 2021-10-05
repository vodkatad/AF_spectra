#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
infile_msi <- args[1]
infile_mss <- args[2]
outfile <- args[3]
colors <- args[4]
image <- args[5]

dmsi <- read.table(infile_msi, sep="\t", quote="", header=TRUE)
dmss <- read.table(infile_mss, sep="\t", quote="", header=TRUE)
dmsi$class <- 'MSI'
dmss$class <- 'MSS'

cbPalette2 <- unlist(strsplit(colors, ','))
d <- rbind(dmsi, dmss)
colnames(d) <- c('mut','estimate','upper','lower','class')
x <- seq(1, nrow(d))
d$x <- as.factor(x)
labels <- substr(d$mut, 2, nchar(as.character(d$mut)))
ggplot(d, aes(x=x, y=estimate, color=class)) +  geom_point(stat="identity", size=3) +
    geom_errorbar(aes(ymin=lower, ymax=upper, x=x, width=0.1, color=class), size=0.6)+theme_bw()+ggtitle('dN/dS')+ylab('ML estimate')+
  theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1),axis.title.y=element_text(size=15))+scale_color_manual(values=cbPalette2)+scale_x_discrete(name="Mutation type", breaks=x,labels=labels)
ggsave(outfile)

if (! is.null(image) ) {
  save.image(image)
}