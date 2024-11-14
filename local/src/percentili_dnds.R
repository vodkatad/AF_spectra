library(ggplot2)

d <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMO_BASALE-tmm.tsv.gz'), sep="\t", header=T, row.names=1)


ave <- apply(d, 1, mean)

qq <- quantile(ave, probs=seq(0,1, 0.1))
genes <- data.frame(gene=names(ave), perc=cut(ave, qq, labels=names(qq)[-1], include.lowest=TRUE), value=ave)


d1 <- read.table('~/mut_CRC0441.tsv', sep="\t", header=T)
d2 <- read.table('~/mut_CRC1430.tsv', sep="\t", header=T)
intersect(d1$gene, d2$gene)
genes$ismut <- as.factor(ifelse(genes$gene %in% d1$gene, 'CRC0441', ifelse(genes$gene %in% d2$gene, 'CRC1430', 'N')))
genes$percc <- as.numeric(gsub('%', '', genes$perc))

ggplot(data=genes, aes(x=ismut, y=percc))+geom_boxplot(outlier.shape=NULL)+geom_jitter(height=0)
ggplot(data=genes, aes(x=ismut, y=value))+geom_boxplot(outlier.shape=NULL)+geom_jitter(height=0)+scale_y_log10()
