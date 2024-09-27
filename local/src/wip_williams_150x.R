library(ggplot2)
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_5.Rdata')

palette_df <- readRDS('/scratch/trcanmed/AF_spectra/local/share/data/palette.rds')
pal <- palette_df$palette
names(pal) <- palette_df$model

our <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/williams_bin_0.12_0.24.r2.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)

our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)

# TODO divide by length - study williams only on CN 1-2-3?
# 1-b/d

bd <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/b_bd.txt', sep="\t", stringsAsFactors = FALSE, header=TRUE)
len <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/len_cn123.txt', sep="\t", stringsAsFactors = FALSE, header=TRUE)
colnames(bd)[1] <- 'sample'

m <- merge(our, len, by="sample")

bd$model <- sapply(bd$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
bd$clone <- sapply(bd$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
bd$model_clone <- paste0(bd$model, "_", bd$clone)


m2 <- merge(m, bd, by="model_clone")

#m2$MRw <- (m2$intercept / m2$len)*(1-m2$d.b)
m2$MRw <- (m2$intercept / m2$len)*(m2$X.b.d..b)
m2$MRw <-  m2$MRw / 0.000000001
  
y_breaks <- guess_ticks(m2$MRw)

ggplot() + 
  geom_point(data=m2, aes(x=model.y, y=MRw, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.7))+
  ylab('Z')+xlab('PDTs')+
  scale_color_manual(values=pal)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   


textSize <- 15
largerSize <- 20
slidetheme <- theme(
  text = element_text(size = textSize, family='sans'),
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(color = "black"),
  panel.background = element_blank()
)

ggplot() + 
  geom_point(data=m2, aes(x=model.x, y=MRw, color=model_clone), stat="identity", position=position_dodge(0.7))+
  ylab('MRw')+xlab('PDTs')+
  scale_color_manual(values=pal)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  slidetheme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   



# which length?

mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_SNV', sep="\t", header=FALSE)
colnames(mredu) <- c('sample', 'MRedu')
# MR_edu
mredu$MRedu <- mredu$MRedu / 0.000000001
mm <- merge(mredu, m2, by="sample", by.y="sample.x")
ggplot() + 
  geom_point(data=mm, aes(x=MRedu, y=MRw, color=model_clone))+
  ylab('MR w')+xlab('MR edu10-9')+
  scale_color_manual(values=pal)+
  slidetheme+theme(legend.position="none",
                   legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   



mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_conte_SNV', sep="\t", header=FALSE)
colnames(mredu) <- c('sample', 'MRedu')
# MR_conte_SNV
mredu$MRedu <- mredu$MRedu / 0.000000001
mm <- merge(mredu, m2, by="sample", by.y="sample.x")
ggplot() + 
  geom_point(data=mm, aes(x=MRedu, y=MRw, color=model_clone))+
  ylab('MRw')+xlab('MR conte10-9')+
  scale_color_manual(values=pal)+
  slidetheme+theme(legend.position="none", 
                   legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   

## aggiungiamo T0
mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_SNV', sep="\t", header=FALSE)
colnames(mredu) <- c('sample', 'MRedu')
m3 <- m2
m3[m3$sample.x=="CRC0282-05-0", 'sample.x'] <-'CRC0282-05-1-C'
m3[m3$sample.x=="CRC1307-09-0", 'sample.x'] <-'CRC1307-09-1-B'
m3[m3$sample.x=="CRC1502-03-0", 'sample.x'] <-'CRC1502-03-1-A'
m3$shape <- ifelse(is.na(m3$clone2), 'T0', 'T1')
# MR_edu
mredu$MRedu <- mredu$MRedu / 0.000000001
mm <- merge(mredu, m3, by="sample", by.y="sample.x")
ggplot() + 
  geom_point(data=mm, aes(x=MRedu, y=MRw, color=model_clone, shape=shape))+
  ylab('MR w')+xlab('MR edu10-9')+
  scale_color_manual(values=pal)+
  slidetheme+theme(legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   



### 150x MR
mr150 <- read.table('/scratch/trcanmed/AF_spectra/dataset150x/MR_edu_SNV', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mr150) <- c('sample', 'MR150edu')
mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_SNV', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mredu) <- c('sample', 'MRedu')
# MR_conte_SNV
mredu$MRedu <- mredu$MRedu / 0.000000001
mr150$MR150edu <- mr150$MR150edu / 0.000000001
mm <- merge(mredu, mr150, by="sample")
mm$model <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
mm$clone <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
mm$clone2 <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
mm$model_clone <- paste0(mm$model, "_", mm$clone)
ggplot() + 
  geom_point(data=mm, aes(x=MRedu, y=MR150edu, color=model_clone))+
  ylab('MR 150x')+xlab('MR 30x')+
  scale_color_manual(values=pal)+
  slidetheme+theme(legend.position="none",
                   legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+ylim(0,10)+xlim(0,10)+
  geom_abline(intercept=0, slope=1)

### 150 univmutect
mr150 <- read.table('/scratch/trcanmed/AF_spectra/dataset150x/MR_edu_univMutect', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mr150) <- c('sample', 'MR150edu')
mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_univMutect', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mredu) <- c('sample', 'MRedu')
# MR_conte_SNV
mredu$MRedu <- mredu$MRedu / 0.000000001
mr150$MR150edu <- mr150$MR150edu / 0.000000001
mm <- merge(mredu, mr150, by="sample")
mm$model <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
mm$clone <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
mm$clone2 <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
mm$model_clone <- paste0(mm$model, "_", mm$clone)
ggplot() + 
  geom_point(data=mm, aes(x=MRedu, y=MR150edu, color=model_clone))+
  ylab('MR clonal 150x')+xlab('MR clonal 30x')+
  scale_color_manual(values=pal)+
  slidetheme+theme(legend.position="none",
                   legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+ylim(0,10)+xlim(0,10)+
  geom_abline(intercept=0, slope=1)

### signatures

s150 <- read.table('/scratch/trcanmed/AF_spectra/dataset150x/vitrovivobulk_heatmap_merged_cosmic.tsv', sep="\t", stringsAsFactors = F, header=T)
s30 <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/vitrovivobulk_heatmap_merged_cosmic.tsv', sep="\t", stringsAsFactors = F, header=T)
s30 <- s30[rownames(s150),]
rownames(s150) <- paste0(rownames(s150), '-150x')
rownames(s30) <- paste0(rownames(s30), '-30x')
ss <- rbind(s150, s30)

setwd('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables')
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_2a_cosmic.pdf.Rdata')
library(RColorBrewer)
library(pheatmap)
palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model
data <- ss
data$model <- unlist(lapply(strsplit(rownames(data),'-'), function(x){ x[1] }))
data$model <- as.character(unlist(lapply(strsplit(data$model,'_'), function(x){ x[1] })))
data$r <- rownames(data)
orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=FALSE)
orderdf$n <- seq(1, nrow(orderdf))
m <- merge(data, orderdf, by="model")
m$seqm <- unlist(lapply(strsplit(m$r,'_'), function(x){ x[2] }))
m <- m[order(m$n, m$seqm),]


data <- m[, colnames(data)]
data$model <- NULL
rownames(data) <- m$r
data$r <- NULL

annot_rows <- data.frame(row.names=rownames(data))

#annot_rows$sample <-  as.factor(unlist(lapply(strsplit(rownames(annot_rows),'_'), function(x){ x[length(x)] })))
annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] }))
annot_rows$model <- as.factor(unlist(lapply(strsplit(annot_rows$model,'_'), function(x){ x[1] })))
# add back LMX P
if (addlm!= "no") {
  annot_rows$model <- paste0(annot_rows$model, ifelse(!grepl('\\d$', annot_rows$model), '', ifelse(annot_rows$model=="CRC0282", 'PR', 'LM')))
  names(pal) <- paste0(names(pal), ifelse(!grepl('\\d$', names(pal)), '', ifelse(names(pal)=="CRC0282", 'PR', 'LM')))
}
colnames(data) <- gsub('X', 'SBS', colnames(data))
data <- t(data)  

annot_colors <- list(model=pal)
ph <- pheatmap(data, show_colnames = TRUE, show_rownames = TRUE,  
               cluster_cols=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  
               color=brewer.pal(9, 'Blues'), breaks=seq(0, 0.4, by=0.05),
               cluster_rows=FALSE, annotation_legend=FALSE)
ph$gtable[[1]][[1]]$children[[1]]$gp$lwd <- 0.001
ph

pheatmap(data, show_colnames = TRUE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  
         color=brewer.pal(5, 'Blues'),
         cluster_rows=FALSE, annotation_legend=FALSE)


for (i in seq(1:nrow(s150))) {
  sa150 <- rownames(s150)[i]
  gg <- strsplit(sa150, '-')[[1]][1]
  sa30 <- rownames(s30)[grepl(gg, rownames(s30))]
  x150v <- as.numeric(s150[sa150,])
  y30v <- as.numeric(s30[sa30,])
  dp <- data.frame(d150x=x150v, d30x=y30v)
  print(ggplot(data=dp, aes(x=d30x, y=d150x))+geom_point()+ggtitle(gg)+slidetheme)
  sp <- cor.test(x150v, y30v, method='spearman')
  print(gg)
  print(c(sp$estimate, sp$p.value))
}

### same clones  HERE HERE WIP SBS 150x ######################

s150 <- read.table('/scratch/trcanmed/AF_spectra/dataset150x/vitrovivobulk_heatmap_merged_cosmic.tsv', sep="\t", stringsAsFactors = F, header=T)
s30 <- read.table('/scratch/trcanmed/AF_spectra/dataset150x/manual_cosmic.tsv', sep="\t", stringsAsFactors = F, header=T)

rownames(s150)[1:3] <- paste0(rownames(s150)[1:3], '-150x')
rownames(s30) <- paste0(rownames(s30), '-30x')
ss <- rbind(s150, s30)

# remove 282
ss <- ss[!grepl('282', rownames(ss)),]
colnames(ss) <- gsub('X', 'SBS', colnames(ss))
pd <- data.frame(SBS1=ss$SBS1, SBS8=ss$SBS8, SBS18=ss$SBS18)
pd$model <- substr(rownames(ss), 0, 7)
pd$class <- ifelse(grepl('bulk', rownames(ss)), 'Pre-existing', ifelse(grepl('30x', rownames(ss)), 'MA 30x', 'MA 150x'))

library(reshape)
pdlong <- melt(data=pd, id=c("class", 'model'))
colnames(pdlong) <- c('class', 'model', 'Signature', 'exp')
pdlong$class <- factor(pdlong$class, c('Pre-existing', 'MA 30x', 'MA 150x'))

ggplot(data=pdlong, aes(y=exp,fill=Signature, x=class))+
  geom_col(color='black')+
  facet_wrap(~model+Signature)+ylab('Relative contribution')+#,expand = c(-0.1, -0.1)
  xlab('')+unmute_theme+scale_y_continuous(breaks=(seq(0, 0.6, by=0.15)),limits=c(-0.05, 0.6),expand = c(0, 0))#+theme(legend.position="bottom",axis.ticks.x = element_blank(),strip.background = element_blank(),strip.text = element_blank())



ggplot(data=pdlong, aes(y=exp,fill=Signature, x=class))+
  geom_col(color='black')+
  facet_wrap(~model+Signature)+ylab('Relative contribution')+#,expand = c(-0.1, -0.1)
  xlab('')+scale_y_continuous(breaks=(seq(0, 0.6, by=0.15)),limits=c(-0.05, 0.6),expand = c(0, 0))+theme_bw(base_size=10)


cor.test(unlist(ss['CRC1307_vitroMA-30x',]), unlist(ss['CRC1307_vitroMA-150x',]), method="spearman")
cor.test(unlist(ss['CRC1502_vitroMA-30x',]), unlist(ss['CRC1502_vitroMA-150x',]), method="spearman")
###############################


setwd('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables')
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_2a_cosmic.pdf.Rdata')
library(RColorBrewer)
library(pheatmap)
palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model
data <- ss
data$model <- unlist(lapply(strsplit(rownames(data),'-'), function(x){ x[1] }))
data$model <- as.character(unlist(lapply(strsplit(data$model,'_'), function(x){ x[1] })))
data$r <- rownames(data)
orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=FALSE)
orderdf$n <- seq(1, nrow(orderdf))
m <- merge(data, orderdf, by="model")
m$seqm <- unlist(lapply(strsplit(m$r,'_'), function(x){ x[2] }))
m <- m[order(m$n, m$seqm),]


data <- m[, colnames(data)]
data$model <- NULL
rownames(data) <- m$r
data$r <- NULL

annot_rows <- data.frame(row.names=rownames(data))

#annot_rows$sample <-  as.factor(unlist(lapply(strsplit(rownames(annot_rows),'_'), function(x){ x[length(x)] })))
annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] }))
annot_rows$model <- as.factor(unlist(lapply(strsplit(annot_rows$model,'_'), function(x){ x[1] })))
# add back LMX P
if (addlm!= "no") {
  annot_rows$model <- paste0(annot_rows$model, ifelse(!grepl('\\d$', annot_rows$model), '', ifelse(annot_rows$model=="CRC0282", 'PR', 'LM')))
  names(pal) <- paste0(names(pal), ifelse(!grepl('\\d$', names(pal)), '', ifelse(names(pal)=="CRC0282", 'PR', 'LM')))
}
colnames(data) <- gsub('X', 'SBS', colnames(data))
data <- t(data)  

annot_colors <- list(model=pal)
ph <- pheatmap(data, show_colnames = TRUE, show_rownames = TRUE,  
               cluster_cols=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  
               color=brewer.pal(9, 'Blues'), breaks=seq(0, 0.4, by=0.05),
               cluster_rows=FALSE, annotation_legend=FALSE)
ph$gtable[[1]][[1]]$children[[1]]$gp$lwd <- 0.001
ph

pheatmap(data, show_colnames = TRUE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  
         color=brewer.pal(5, 'Blues'),
         cluster_rows=FALSE, annotation_legend=FALSE)


for (i in seq(1:nrow(s150))) {
  sa150 <- rownames(s150)[i]
  gg <- strsplit(sa150, '-')[[1]][1]
  sa30 <- rownames(s30)[grepl(gg, rownames(s30))]
  x150v <- as.numeric(s150[sa150,])
  y30v <- as.numeric(s30[sa30,])
  dp <- data.frame(d150x=x150v, d30x=y30v)
  print(ggplot(data=dp, aes(x=d30x, y=d150x))+geom_point()+ggtitle(gg)+slidetheme)
  sp <- cor.test(x150v, y30v, method='spearman')
  print(gg)
  print(c(sp$estimate, sp$p.value))
}


### 150x MR 60x
mr150 <- read.table('/scratch/trcanmed/AF_spectra/dataset150x_60x/MR_edu_SNV', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mr150) <- c('sample', 'MR150edu')
mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_SNV', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mredu) <- c('sample', 'MRedu')
# MR_conte_SNV
mredu$MRedu <- mredu$MRedu / 0.000000001
mr150$MR150edu <- mr150$MR150edu / 0.000000001
mm <- merge(mredu, mr150, by="sample")
mm$model <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
mm$clone <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
mm$clone2 <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
mm$model_clone <- paste0(mm$model, "_", mm$clone)
ggplot() + 
  geom_point(data=mm, aes(x=MRedu, y=MR150edu, color=model_clone))+
  ylab('MR 150x')+xlab('MR 30x')+
  scale_color_manual(values=pal)+
  slidetheme+theme(legend.position="none",
                   legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+ylim(0,10)+xlim(0,10)+
  geom_abline(intercept=0, slope=1)

## seghe
mr150 <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_binSNV', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mr150) <- c('sample', 'MR150edu')
mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_SNV', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mredu) <- c('sample', 'MRedu')
# MR_conte_SNV
mredu$MRedu <- mredu$MRedu / 0.000000001
mr150$MR150edu <- mr150$MR150edu / 0.000000001
mm <- merge(mredu, mr150, by="sample")
mm$model <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
mm$clone <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
mm$clone2 <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
mm$model_clone <- paste0(mm$model, "_", mm$clone)
ggplot() + 
  geom_point(data=mm, aes(x=MRedu, y=MR150edu, color=model_clone))+
  ylab('MR 150x')+xlab('MR 30x')+
  scale_color_manual(values=pal)+
  slidetheme+theme(legend.position="none",
                   legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+ylim(0,10)+xlim(0,10)+
  geom_abline(intercept=0, slope=1)



mr150 <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_univMutect', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mr150) <- c('sample', 'MR150edu')
mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_SNV', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(mredu) <- c('sample', 'MRedu')
# MR_conte_SNV
mredu$MRedu <- mredu$MRedu / 0.000000001
mr150$MR150edu <- mr150$MR150edu / 0.000000001
mm <- merge(mredu, mr150, by="sample")
mm$model <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
mm$clone <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
mm$clone2 <- sapply(mm$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
mm$model_clone <- paste0(mm$model, "_", mm$clone)
ggplot() + 
  geom_point(data=mm, aes(x=MRedu, y=MR150edu, color=model_clone))+
  ylab('MR 150x')+xlab('MR 30x')+
  scale_color_manual(values=pal)+
  slidetheme+theme(legend.position="none",
                   legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+ylim(0,10)+xlim(0,10)+
  geom_abline(intercept=0, slope=1)



