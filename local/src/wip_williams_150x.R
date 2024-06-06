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

m2$MRw <- (m2$intercept / m2$len)*m2$d.b
m2$MRw <-  m2$MRw / 0.000000001
  
y_breaks <- guess_ticks(m2$MRw)

ggplot() + 
  geom_point(data=m2, aes(x=model, y=MRw, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.7))+
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
# MR_conte_SNV
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
