d <- read.table("/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/signatures/all_cosmic_fit.tsv", sep="\t", header=T, row.names=1)

models <- unique(substr(colnames(d), 0, 7))

get_sign <- function(model, data, sign) {
  sub <- data[sign, grepl(model, colnames(data))]
  return(c(sub[, grepl('PRX', colnames(sub))], sub[, grepl('LMX', colnames(sub))]))
}

sbs1 <- as.data.frame(t(sapply(models, get_sign, d, 'Signature.1')))
sbs8 <- as.data.frame(t(sapply(models, get_sign, d, 'Signature.8')))
sbs18 <- as.data.frame(t(sapply(models, get_sign, d, 'Signature.18')))
colnames(sbs1) <- c('PRX', 'LMX')
colnames(sbs8) <- c('PRX', 'LMX')
colnames(sbs18) <- c('PRX', 'LMX')

library(ggplot2)


pairs <- function(da, title) {
  da$id <- row.names(da)
  m <- melt(da)
  pd <- position_dodge(width=0.2) # inv
  print(ggplot(data=m, aes(x=variable, y=value, group=id))+geom_jitter(position=pd)+geom_line(aes(group=id), position=pd)+ggtitle(title)+theme(axis.text=element_text(size=15)))
}
pairs(sbs1, 'SBS1')
pairs(sbs8, 'SBS8')
pairs(sbs18, 'SBS18')


d$id <- NULL
d <- d[seq(1, nrow(d)-1),]
pheatmap(d)

annot <- data.frame(row.names=colnames(d), class=ifelse(grepl('PRX', colnames(d)), 'PRX', 'LMX'))
pheatmap(d, annotation_col=annot, show_colnames = F)


d2 <- d[rownames(d) %in% c('Signature.8','Signature.18'),]
d2 <-d2[,order(colnames(d2))]
pheatmap(d2, annotation_col=annot, show_colnames = F, cluster_cols = F)
