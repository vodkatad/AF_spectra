MR_f  <- snakemake@input[['MR']]
MRuniv_f  <- snakemake@input[['MRclo']]

colors <- snakemake@input[['palette']]
model_colors <- snakemake@input[['palette2']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

palette_df <- readRDS(model_colors)
pal2 <- palette_df$palette
names(pal2) <- palette_df$model

subclonal <- read.table(MRuniv_f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(subclonal) <- c('sample', 'MR')
subclonal$MR_edu_cl <- subclonal$MR / 0.000000001
MR <- read.table(MR_f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(MR) <- c('sample', 'MR')
MR$MR_edu <- MR$MR / 0.000000001

# remove MSI
MR <- MR[!grepl('CRC0282', MR$sample),]
subclonal <- subclonal[!grepl('CRC0282', subclonal$sample),]

merged <- merge(subclonal, MR, by='sample')
merged$smodel <- substr(merged$sample, 0, 7)
#merged <- merged[merged$smodel != 'CRC0282',]

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
#p <- ggplot() +
#  geom_point(data=merged, aes(x=MR_edu, y=MR_edu_cl, color=model_clone), stat="identity", size=2)+
#  geom_point(data=pdata, aes(x=mean.x, y=mean.y), size=2)+
#  geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, ymin=lower.y, ymax=upper.y), width=.05, size=.2)+
#  geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, xmin=lower.x, xmax=upper.x), width=.05, size=.2)+
#  scale_fill_manual(values=pal2)+scale_color_manual(values=pal, guide='none')+
  #xlab('# SNVs gained')+ylab('# clonal SNVs gained')+theme(legend.position="right")+geom_abline(slope=1, intercept=0)
#  xlab('MR')+ylab('MR with clonal calls')+theme(legend.position="right")+geom_abline(slope=1, intercept=0)+
#  unmute_theme
xb <- guess_ticks(values=merged$MR_edu, fixed_max=4)
yb <- guess_ticks(values=merged$MR_edu_cl, fixed_max=4)


p <- ggplot() +
  geom_point(data=merged, aes(x=MR_edu, y=MR_edu_cl, color=model_clone),  stat="identity", size=2, shape=16, alpha=0.8)+
  scale_color_manual(values=pal, guide="none")+
  scale_y_continuous(breaks=yb, limits=c(0, max(yb)), expand = c(0, 0))+
  scale_x_continuous(breaks=xb, limits=c(0, max(xb)), expand = c(0, 0))+
  #geom_point(data=pdata, aes(x=mean.x, y=mean.y, fill=model), color="black", alpha=0.5, size=1.5, shape=22, stroke=0.2)+ scale_fill_manual(values=pal2, guide="none")+
  #geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, ymin=lower.y, ymax=upper.y), width=.05, size=.2)+
  #geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, xmin=lower.x, xmax=upper.x), width=.05, size=.2)+
  xlab('MR')+ylab('MR with clonal calls')+theme(legend.position="right")+geom_abline(slope=1, intercept=0)+
  unmute_theme


death_conversion_dpi96 <- 96/72
ggsave(outplot, plot=p, width=60*death_conversion_dpi96, height=60*death_conversion_dpi96, units="mm")

sink(log_f)
cor.test(pdata$mean.x, pdata$mean.y, method='spearman')
print('all')
p <- cor.test(merged$MR_edu, merged$MR_edu_cl, method='spearman')
print(p)
print(p$p.value)
sink()

save.image(paste0(outplot, '.Rdata'))
