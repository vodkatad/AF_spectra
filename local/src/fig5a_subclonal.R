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
  geom_point(aes(color=model_clone, shape=PDT, fill=model_clone), stat="identity", size=1, position=position_dodge(0.2))+
  unmute_theme+scale_color_manual(values=pal, guide="none")+scale_shape_manual(values=c(24,25,22,23))+#scale_shape_manual(values=c(18,23,20,19))+
  scale_fill_manual(values=pal, guide="none")+
  xlab('MR')+ylab('# subclonal SNVs')+theme(legend.position="right")

# average errors ############
confidence_interval <- function(vector, interval) {
   # Standard deviation of sample
   vec_sd <- sd(vector)
   # Sample size
   n <- length(vector)
   # Mean of sample
   vec_mean <- mean(vector)
   # Error according to t distribution
   error <- qt((1+interval)/2, df = n - 1) * vec_sd / sqrt(n)
   # Confidence interval as a vector
   result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
   return(result)
}
 
LEVEL <- 0.99
ic_clones <- sapply(unique(merged$model), function(x) { confidence_interval(merged[merged$model==x,'MR_edu'], LEVEL) })
colnames(ic_clones) <- unique(merged$model)
pdataMR <- as.data.frame(t(ic_clones))
pdataMR$model <- rownames(pdataMR)


ic_clones <- sapply(unique(merged$model), function(x) { confidence_interval(merged[merged$model==x,'n'], LEVEL) })
colnames(ic_clones) <- unique(merged$model)
pdatasub <- as.data.frame(t(ic_clones))
pdatasub$model <- rownames(pdatasub)

pdata <- merge(pdataMR, pdatasub, by='model')
# 
p <- ggplot() +
   geom_point(data=merged, aes(x=MR_edu, y=n, color=model_clone, shape=PDT, fill=model_clone), stat="identity", size=1, position=position_dodge(0.2))+
   unmute_theme+scale_color_manual(values=pal, guide="none")+scale_shape_manual(values=c(24,25,22,23))+#scale_shape_manual(values=c(18,23,20,19))+
   geom_point(data=pdata, aes(x=mean.x, y=mean.y), shape=1, size=2)+
   geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, ymin=lower.y, ymax=upper.y), width=.05, size=.2)+
   geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, xmin=lower.x, xmax=upper.x), width=.05, size=.2)+
   scale_fill_manual(values=pal, guide="none")+
   xlab('MR')+ylab('# subclonal SNVs')+theme(legend.position="right")


# pearson pvalue
sink(log_f)
for (me in c('spearman', 'pearson')) {
  print('all')
  print(cor.test(merged$MR_edu, merged$n, method=me))
  print(nrow(merged))
  for (m in unique(merged$model)) {
    print(m)
    subm <- merged[merged$model == m, ]
    print(cor.test(subm$MR_edu, subm$n, method=me))
    print(nrow(subm))
  }
  print('avg')
  print(cor.test(pdata$mean.x, pdata$mean.y, method=me))
  print(nrow(pdata))
}
sink()


ggsave(outplot, plot=p, width=89, height=56, units="mm")
save.image(paste0(outplot, '.Rdata'))