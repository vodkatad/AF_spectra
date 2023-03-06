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
output_mrca <- args[5]
mr_f <- args[6]

#files <- list.files(path="./", pattern = paste0("*_",wanted,".ovcnokdelta.tsv.gz"), recursive = F, full.names = F)
files <- c()
for (i in seq(7, length(args))) {
  files <- c(files, args[i])
}

if (length(files)==1) {
  file.create(output_plot) 
  file.create(output_muts) # fixme generate files with 0 or .. ?
  file.create(output_n)
  private_story <- data.frame(matrix(NA, nrow=1, ncol=5))
  colnames(private_story) <- c('clone_ref', 'clone_other', 'private_gen', 'private_mut_ref', 'common_mut')
  write.table(private_story, file=output_mrca, quote=F, sep="\t")
  q(save="no", status=0)
}

load_gain <- function(filename) {
  d <- read.table(filename, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  d <- d[d$V16=="gain",]
  muts <- d$V1
  if (any(grepl('@', muts))) {
    muts <- sapply(strsplit( muts, '@'), function(x) {x[[1]]})
  }
  return(muts)
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

#> table(counts$Freq)

#   1    2    3 
#5274 2091  393 

save.image('pippo.Rdata')
merged <- lapply(gains, function(x) { m <- counts[counts$ugains %in% x,]; return(m)})
merged_named <- lapply(names(merged), function(x) { y <- merged[[x]]; counts <- as.data.frame(table(y$Freq)); counts$name <- x; return(counts) })
pdata <- do.call(rbind, merged_named)
pdata$name <- unlist(lapply(strsplit(pdata$name,'_'), function(x) {x[[1]][1]}))
colnames(pdata) <- c("shared_by","count","sample")
#max(as.numeric(as.character(pdata$shared_by)))
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

## load MR and process end clones in pairs
mr <- read.table(mr_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
end_clones <- names(gains)
private_story <- data.frame(matrix(NA, nrow=length(end_clones)**2-length(end_clones), ncol=5))
colnames(private_story) <- c('clone_ref', 'clone_other', 'private_gen', 'private_mut_ref', 'common_mut')
# we do both orders cause we compute unique muts of each clone separately, i is the one we refer to (and get MR of)
k <- 1
for (i in seq(1, length(end_clones))) {
  reference_clone <- strsplit(end_clones[i], "_")[[1]][1]
  for (j in seq(1, length(end_clones))) {
    if (i != j) {
      other_clone <- strsplit(end_clones[j], "_")[[1]][1]
      ref_mr <- mr[mr$end == reference_clone, ]
      ref_only_gained <- length(setdiff(gains[[end_clones[i]]], gains[[end_clones[j]]]))
      common_gained <- length(intersect(gains[[end_clones[i]]], gains[[end_clones[j]]]))
      private_story[k,] <- c(reference_clone, other_clone, ref_only_gained/(ref_mr[,'MR_EDU']*ref_mr[,'len_cnok']), ref_only_gained, common_gained)
      k <- k + 1
    }
  }
}
write.table(private_story, file=output_mrca, quote=F, sep="\t")