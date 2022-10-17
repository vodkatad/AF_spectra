#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
outfile <- args[1]
colores <- args[2]



# per ogni input faccio media dopo aver rimosso gli invivo
summarized <- data.frame(model=c(), EDU=c(), cellcounts=c())
for (j in seq(2, length(args))) {
  data <- read.table(args[j], sep="\t", header=TRUE, stringsAsFactors=FALSE)
  data <- data[!grepl('-M', data$end),]
  summarized[j-1, 'model'] <- unique(sapply(data$end, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])}))
  summarized[j-1, 'EDU'] <- mean(data$MR_EDU) / 0.000000001
  summarized[j-1, 'cellcounts'] <- mean(data$MR_conte) / 0.000000001
}
print(summarized)

ggplot(summarized, aes(x=model, color=model))+
geom_errorbar(aes(ymin=EDU, ymax=cellcounts), width=.2, size=1, position=position_dodge(1))+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_color_manual(values=cbPalette)

ggsave(outfile)