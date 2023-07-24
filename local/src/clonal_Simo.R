setwd('/scratch/trcanmed/AF_spectra/datasetV2')

palette_df <- readRDS("/scratch/trcanmed/AF_spectra/local/share/data/palette.rds")
pal <- palette_df$palette
names(pal) <- palette_df$model_clone

palette_df <- readRDS("/scratch/trcanmed/AF_spectra/local/share/data/model_palette.rds")
pal2 <- palette_df$palette
names(pal2) <- palette_df$model_clone

subclonal <- read.table('vitro_gained_binSNV', header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(subclonal) <- c('sample', 'MR')
#subclonal$MR_edu_cl <- subclonal$MR / 0.000000001
subclonal$MR_edu_cl <- subclonal$MR 
MR <- read.table('vitro_gained_SNV', header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(MR) <- c('sample', 'MR')
#MR$MR_edu <- MR$MR / 0.000000001
MR$MR_edu <- MR$MR 

merged <- merge(subclonal, MR, by='sample')
merged$smodel <- substr(merged$sample, 0, 7)
#wanted <- c('CRC1599', 'CRC1307', 'CRC1078')
merged <- merged[merged$smodel != 'CRC0282',]

merged$model <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
merged$PDT <- paste0(merged$model, ifelse(grepl('\\d$', merged$model), ifelse(grepl('282', merged$model), 'PR', 'LM'), ''))
#merged$PDT <- factor(merged$PDT, levels=c('CRC1078LM','CRC1307LM','CRC1599PR','CRC1599LM'))

merged$clone <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
merged$clone2 <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
merged$model_clone <- paste0(merged$model, "_", merged$clone)
merged$time <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

# p <- ggplot(data=merged, aes(x=MR_edu, y=n)) +
#   geom_point(aes(color=model_clone, fill=model_clone), stat="identity", size=1, position=position_dodge(0.2))+
#   unmute_theme+scale_color_manual(values=pal, guide="none")+
#   scale_fill_manual(values=pal, guide="none")+
#   xlab('MR')+ylab('# subclonal SNVs')+theme(legend.position="right")

# average errors ############
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- sd(vector)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

LEVEL <- 0.99
ic_clones <- sapply(unique(merged$model), function(x) { confidence_interval(merged[merged$model==x,'MR_edu'], LEVEL) })
colnames(ic_clones) <- unique(merged$model)
pdataMR <- as.data.frame(t(ic_clones))
pdataMR$model <- rownames(pdataMR)


ic_clones <- sapply(unique(merged$model), function(x) { confidence_interval(merged[merged$model==x,'MR_edu_cl'], LEVEL) })
colnames(ic_clones) <- unique(merged$model)
pdatasub <- as.data.frame(t(ic_clones))
pdatasub$model <- rownames(pdatasub)

pdata <- merge(pdataMR, pdatasub, by='model')
# 
ggplot() +
  geom_point(data=merged, aes(x=MR_edu, y=MR_edu_cl, color=model_clone), stat="identity", size=2)+
  theme_bw()+
  geom_point(data=pdata, aes(x=mean.x, y=mean.y, shape=model, fill=model), size=2)+
  geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, ymin=lower.y, ymax=upper.y), width=.05, size=.2)+
  geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, xmin=lower.x, xmax=upper.x), width=.05, size=.2)+
  scale_fill_manual(values=pal2)+scale_color_manual(values=pal, guide='none')+
  xlab('# SNVs gained')+ylab('# clonal SNVs gained')+theme(legend.position="right")
