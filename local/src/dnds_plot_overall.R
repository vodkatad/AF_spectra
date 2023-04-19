#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
colors <- args[3]
image <- args[4]

d <- read.table(infile, sep="\t", quote="", header=FALSE)
colnames(d) <- c('model','estimate', 'upper', 'lower')
d <- d[d$model != "overall",]

cbPalette2 <- unlist(strsplit(colors, ','))


ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)


ggplot(d, aes(x=model, y=estimate, color=model)) +  geom_point(stat="identity", size=2) +
    geom_errorbar(aes(ymin=lower, ymax=upper, x=model, width=0.1, color=model), size=0.2)+theme_bw()+ggtitle('dN/dS nonsyn+truncating')+ylab('ML estimate')+xlab('')+
  ctheme+scale_color_manual(values=cbPalette2)

ggsave(outfile)
if (!is.na(image) && image != "") {
  save.image(image)
}

q(save="no")

#load('dndsvitro_overall.png.Rimage')
d <- d[!grepl('2nd', d$model),]
cbPalette2 <- unique(cbPalette2)
ggplot(d, aes(x=model, y=estimate, color=model)) +  geom_point(stat="identity", size=3) +
    geom_errorbar(aes(ymin=lower, ymax=upper, x=model, width=0.1, color=model), size=0.6)+theme_bw()+ggtitle('dN/dS nonsyn+truncating')+ylab('ML estimate')+xlab('')+
  ctheme+scale_color_manual(values=cbPalette2)
ggsave('dndsvitro_overall_simpler.png')
#savehistory()