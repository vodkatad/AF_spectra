#!/usr/bin/env Rscript
library(RColorBrewer)
library(ggplot2)
library(scales)
#setwd("/mnt/trcanmed/snaketree/prj/AF_spectra/dataset/CRC1307_platypus_nobin/")
args <- commandArgs(trailingOnly = T)
wanted <- args[1]
output_plot <- args[2]
output_n <- args[3]
output_muts <- args[4]


files <- list.files(path="./", pattern = paste0("*_",wanted,".ovcnokdelta.tsv.gz"), recursive = F, full.names = F)

load_gain <- function(filename) {
  d <- read.table(filename, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  d <- d[d$V16=="gain",]
  return(d[,1])
}

gains <- sapply(files, load_gain)

ugains <- unlist(gains)
counts <- as.data.frame(table(ugains))
#counts$class <- ifelse(counts$Freq == 3)
# 
col1 <- rev(tail(brewer.pal(n = 9, name = "Blues"), n=3))
col2 <- rev(tail(brewer.pal(n = 9, name = "Reds"), n=9))
#col3 <- tail(brewer.pal(n = 9, name = "BuPu"), n=1)
col4 <- rev(tail(brewer.pal(n = 9, name = "Greens"), n=3))
#col5 <- brewer.pal(n = 6, name = "Set3")
colo <- c(col1, col4, col2)


merged <- lapply(gains, function(x) { m <- counts[counts$ugains %in% x,]; return(m)})
merged_named <- lapply(names(merged), function(x) { y <- merged[[x]]; counts <- as.data.frame(table(y$Freq)); counts$name <- x; return(counts) })
pdata <- do.call(rbind, merged_named)
pdata$name <- unlist(lapply(strsplit(pdata$name,'_'), function(x) {x[[1]][1]}))
colnames(pdata) <- c("shared_by","count","sample")
#max(as.numeric(as.character(pdata$shared_by)))
save.image('pippo.Rdata')
#pdata$shared_by <- factor(pdata$shared_by, levels = c(1,2,3,4,5,6,7,11)) 
#order of levels need to be guaranteed:
l <- as.numeric(as.character(levels(pdata$shared_by)))
l <- l[order(l)]
pdata$shared_by <- factor(pdata$shared_by, levels=l)
if (length(levels(pdata$shared_by)) < length(colo)) {
  ggplot(pdata,aes(x = sample, y = count, fill = shared_by)) +geom_bar(position = "fill",stat = "identity")+scale_fill_manual(values=colo) +scale_y_continuous(labels = percent_format())+theme_bw()+theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(wanted)
  ggsave(output_plot)
} else {
  ggplot(pdata,aes(x = sample, y = count, fill = shared_by)) +geom_bar(position = "fill",stat = "identity")+scale_fill_manual(values=c(colo, colo)) +scale_y_continuous(labels = percent_format())+theme_bw()+theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(wanted)
  ggsave(output_plot)
}
write.table(pdata, file=output_n, quote=F, sep="\t")


### specific 1-2-3 muts
n_common <- do.call(rbind, merged)
# those found in more than 1 sample are repeated:
u_n_common <- unique(n_common)
write.table(u_n_common, file=output_muts, quote=F, sep="\t")