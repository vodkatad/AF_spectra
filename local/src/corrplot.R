#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
toremove <- args[1]
colors <- args[2]
infile_x <- args[3]
infile_y <- args[4]
outfile <- args[5]
column <- args[6]
name_x <- args[7]
name_y <- args[8]
rate <- args[9]

cbPalette <- unlist(strsplit(colors, ','))
save.image('pluto.Rdata')
dx <- read.table(infile_x, sep="\t", header=TRUE, stringsAsFactors=FALSE)
dy <- read.table(infile_y, sep="\t", header=TRUE, stringsAsFactors=FALSE)
dx <- dx[, colnames(dx) %in% c('model',column), drop=FALSE]
dy <- dy[, colnames(dy) %in% c('model',column), drop=FALSE]
m <- merge(dx, dy, by='model')
colnames(m) <- c('model',name_x, name_y)
m <- m[m$model!=toremove,]
if (rate == "yes") {
  m[,name_y] <- m[,name_y] / 0.000000001
  m[,name_x] <- m[,name_x] / 0.000000001
}
ci <- cor.test(m[, name_x], m[, name_y])

ggplot(m, aes_string(x=name_x, y=name_y)) +  geom_point(aes(color=model), size=3) + geom_smooth(method='lm')+
  theme_bw()+labs(caption=paste0('pearson=', ci$estimate, ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=cbPalette)
#ggplot(m, aes(x=MR_SNV, y=MR_indel)) +  geom_point(aes(color=model), size=3) + geom_smooth(method='lm', se=TRUE)+
#  theme_bw()+labs(caption=paste0('pearson=', ci$estimate, ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=cbPalette)
ggsave(outfile)