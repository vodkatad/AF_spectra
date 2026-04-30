load('/scratch/trcanmed/AF_spectra/datasetV2/p.Rdata')
library(ggpubr)

cos_sim_ori_rec <- cos_sim_matrix(mut_mat, ff$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$s = row.names(cos_sim_ori_rec)
cos_sim_ori_rec$sample <- unlist(lapply(strsplit(cos_sim_ori_rec$s,"_"), function(x) x[length(x)]))
cos_sim_ori_rec$n <- seq(1, nrow(cos_sim_ori_rec))
ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=reorder(s, n), fill=sample)) +
  geom_bar(stat="identity") +
  coord_flip(ylim=c(0.7,1)) +
  ylab("Cosine similarity\n original VS reconstructed") + 
  xlab("") +theme_bw()+ theme(axis.text.y = element_text(size=10), axis.title.x=element_text(size=15))+
  geom_hline(aes(yintercept=.95))+theme(legend.position="none")



data <- as.matrix(ff$contribution)
#select <- which(rowSums(data) > 10)
data = t(data)
data = data/rowSums(data)
#select <- which(colSums(data)>1)
#data <- data[,select]
annot_rows <- data.frame(row.names=rownames(data))

#annot_rows[grepl('M',rownames(annot_rows)), 'sample'] <- 'vivo'
annot_rows$time <-  as.factor(unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[[3]] })))
annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] }))
annot_rows$model <- as.factor(unlist(lapply(strsplit(annot_rows$model,'_'), function(x){ x[1] })))
#annot_cols <- data.frame(signature=c('Clock','Clock','ROS','smoke','MMR','MMR','NA','NA'), row.names=colnames(data))

# cosmic specific XXX WARNING
colnames(data) <- unlist(lapply(strsplit(colnames(data),'.', fixed=TRUE), function(x){ x[length(x)] }))

pheatmap(data, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE,  annotation_row=annot_rows, color=brewer.pal(9, 'PuBu'))

#
data2 <- data[!grepl('-M', rownames(data)),]
annot_rows2 <- annot_rows[annot_rows$time %in% c(0, 1),]
pheatmap(data2, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE,  annotation_row=annot_rows2, color=brewer.pal(9, 'PuBu'))


data3 <- data2[,apply(data2, 2, sum) > 1 ]
pheatmap(data3, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE,  annotation_row=annot_rows2, color=brewer.pal(9, 'PuBu'))

###################
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_10.Rdata')
pali <- readRDS("/scratch/trcanmed/AF_spectra/local/share/data/IANG_allclones_palette.rds")
palma <- readRDS("/scratch/trcanmed/AF_spectra/local/share/data/palette.rds")
palette_df <- rbind(palma, pali)
pal <- palette_df$palette
names(pal) <- palette_df$model
names(pal) <- gsub('-', '_', names(pal))

dd <- data.frame(ratio=data2[,8] / data2[,1], row.names= rownames(data2))
dd$time <-  as.factor(unlist(lapply(strsplit(rownames(dd),'-'), function(x){ x[[3]] })))
dd$model <-  as.factor(unlist(lapply(strsplit(rownames(dd),'-'), function(x){ x[[1]] })))
dd$clone <-  as.factor(unlist(lapply(strsplit(rownames(dd),'-'), function(x){ x[[2]] })))

# only in vivo T1 for this
dd <- dd[rownames(dd) != 'CRC1599LM-07-0',]
pal <- RColorBrewer::brewer.pal(name = 'OrRd', n=8)
pl <- c()
for (m in unique(dd$model)) {
  md <- dd[dd$model == m, ]
  md0 <- md[md$time == 0,]
  md0$lineage <- NA
  md0$lineage <- paste0('clone', md0[,'clone'])
  md1 <- md[md$time == 1,]
  md1$lineage <- NA
  for (c in unique(md1$clone)) {
    ll <- median(md1[md1$clone == c,'ratio'])
    if (!any(ll ==  md1[md1$clone == c,'ratio'])) {
      # solo CRC1599PR
      ll = max(md1[md1$clone == c,'ratio'])
    }
    md1[which(md1$clone==c & md1$ratio == ll),]$lineage <- paste0('clone', c)
  }
  md <- rbind(md0, md1)
  lenc <- length(unique(md$clone))
  pl <- c(pl, ggplot(data=md, aes(x=time, y= ratio))+geom_line(data=md[!is.na(md$lineage),], aes(group=lineage, color=clone))+
        geom_point(aes(color=clone))+ theme_bw(base_size = 15)+
        scale_color_manual(values=pal[seq(8, 8 - lenc)])++ theme(legend.position="none"))
}

pl <- c()
for (m in unique(dd$model)) {
  
  md <- dd[dd$model == m, ]
  md0 <- md[md$time == 0,]
  #md0$lineage <- NA
  #md0$lineage <- paste0('clone', md0[,'clone'])
  md1 <- md[md$time == 1,]
  #md1$lineage <- NA
  md1 <- md1 %>% summarize(ratio=mean(ratio), .by=c(clone))
  md1$time <- 1
  md1$model <- m
  md1 <- md1[, colnames(md0)]
  md <- rbind(md0, md1)
  md$col <- paste0(md$model,'_', md$clone)
  lenc <- length(unique(md$clone))
  mpal <- pal[names(pal) %in% unique(md$col)]
  pl <- c(pl, ggplot(data=md, aes(x=time, y= ratio))+geom_line(aes(group=clone, color=col))+geom_point(aes(color=col))+theme_bw(base_size = 15)+ggtitle(m)+
            scale_color_manual(values=mpal)+ theme(legend.position="none"))
}


ggarrange(plotlist=pl)

## get estimates
# (deltaTumor * nDiv) / deltaMA

gens <- read.table('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/edt3_MR.tsv', sep="\t", header=T)
gens$time <-  as.factor(unlist(lapply(strsplit(gens$id,'-'), function(x){ x[[3]] })))
gens$model <-  as.factor(unlist(lapply(strsplit(gens$id,'-'), function(x){ x[[1]] })))
gens$clone <-  as.factor(unlist(lapply(strsplit(gens$id,'-'), function(x){ x[[2]] })))

gens <- gens[gens$time != 2,]

# deltaMA : deltaTumor =  ngenMA : ngenTumor
# ngenTumor = (deltaTumor * ngenMA ) / deltaMA
# deltaTumor : ngenTumor = deltaMA : ngenMA

models <- unique(dd$model)
ages <- c()
for (m in models) {
  md <- dd[dd$model == m, ]
  md0 <- md[md$time == 0,]
  md1 <- md[md$time == 1,]
  t1s <- c()
  t0s <- c()
  for (c in unique(md$clone)) {
    ave_t1 <- mean(md1[md1$clone==c,'ratio'])
    t0 <- md0[md0$clone==c,'ratio']
    t0s <- c(t0s, t0)
    t1s <- c(t1s, ave_t1)
  }
  n_gen <- mean(unique(gens[gens$model == m, 'Generations_EDU']))
  ages <- c(ages, (n_gen*mean(t0s))/(mean(t1s)-mean(t0s)) )
}

dfage <- data.frame(row.names=models, ages=ages)

### different plot with SBS1 and 8 as separate lines
library(dplyr)

###################
dd2 <- data.frame(SBS8=data2[,8], SBS1=data2[,1], row.names= rownames(data2))
dd2$time <-  as.factor(unlist(lapply(strsplit(rownames(dd2),'-'), function(x){ x[[3]] })))
dd2$model <-  as.factor(unlist(lapply(strsplit(rownames(dd2),'-'), function(x){ x[[1]] })))
dd2$clone <-  as.factor(unlist(lapply(strsplit(rownames(dd2),'-'), function(x){ x[[2]] })))

# only in vivo T1 for this
dd2 <- dd2[rownames(dd2) != 'CRC1599LM-07-0',]

for (m in unique(dd2$model)) {
  md <- dd2[dd2$model == m, ]
  md0 <- md[md$time == 0,]
  md0$lineage <- NA
  md0$lineage <- paste0('clone', md0[,'clone'])
  md1 <- md[md$time == 1,]
  md1$lineage <- NA
  for (c in unique(md1$clone)) {
    ll <- median(md1[md1$clone == c,'ratio'])
    if (!any(ll ==  md1[md1$clone == c,'ratio'])) {
      # solo CRC1599PR
      ll = max(md1[md1$clone == c,'ratio'])
    }
    md1[which(md1$clone==c & md1$ratio == ll),]$lineage <- paste0('clone', c)
  }
  md <- rbind(md0, md1)
  print(ggplot(data=md, aes(x=time, y= ratio))+geom_line(data=md[!is.na(md$lineage),], aes(group=lineage, color=clone))+geom_point(aes(color=clone))+theme_bw(base_size = 15)+ggtitle(m))
}

# fare vs media e togliere legenda

md$model <- NULL
pd <- melt(md, id=c('time', 'clone'), value=c('SBS1', 'SBS2'))
pdt1 <- pd[pd$time == 1,]
pdt0 <- pd[pd$time == 0,]
#pdt1av <- sapply(unique(pdt1$group), function(x) { su <- pdt0[pdt0$group==x, 'value']; mean(su)})
pdt1 <- pdt1 %>% summarize(value=mean(value), .by=c(clone, variable))
pdt1$time <- 1
pdt1 <- pdt1[, colnames(pdt0)]
pd2 <- rbind(pdt0, pdt1)
pd2$group <- paste0(pd2$clone,pd2$variable)
ggplot(data=pd2,aes(x=time, y=value))+geom_point(aes(color=variable))+geom_line(aes(group=variable))

## che misero!
