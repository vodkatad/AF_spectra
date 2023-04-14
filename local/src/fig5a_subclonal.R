subclonal_f  <- snakemake@input[['subclonal']]
MR_f  <- snakemake@input[['MR']]
colors <- snakemake@input[['colors']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))
library(ggplot2)
library(ggpubr)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model_clone

subclonal <- read.table(subclonal_f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(subclonal) <- c('sample', 'n')
MR <- read.table(MR_f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(MR) <- c('sample', 'MR')
MR$MR_edu <- MR$MR / 0.000000001

merged <- merge(subclonal, MR, by='sample')
merged$smodel <- substr(merged$sample, 0, 7)
wanted <- c('CRC1599', 'CRC1307', 'CRC1078')
merged <- merged[merged$smodel %in% wanted,]

merged$model <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
merged$PDT <- paste0(merged$model, ifelse(grepl('\\d$', merged$model), 'LM', ''))
merged$PDT <- factor(merged$PDT, levels=c('CRC1078LM','CRC1307LM','CRC1599PR','CRC1599LM'))

merged$clone <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
merged$clone2 <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
merged$model_clone <- paste0(merged$model, "_", merged$clone)
merged$time <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

p <- ggplot(data=merged, aes(x=MR_edu, y=n)) +
  geom_point(aes(color=model_clone, shape=PDT, fill=model_clone), stat="identity", size=2, position=position_dodge(0.2))+
  unmute_theme+scale_color_manual(values=pal, guide="none")+scale_shape_manual(values=c(18,18,20,20))+#scale_shape_manual(values=c(18,23,20,19))+
  scale_fill_manual(values=pal, guide="none")+
  xlab('MR, mut/(division*bp) *10^-9')+ylab('# subclonal SNVs')


ggsave(outplot, plot=p, width=89, height=56, units="mm")
save.image(paste0(outplot, '.Rdata'))