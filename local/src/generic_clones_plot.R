#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
colors <- args[3]
yname <- args[4]

our <- read.table(infile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(our) <- c('sample', 'MR_edu')
our <- our[!is.na(our$MR_edu),]
#our$MR <- our$MR_edu / 0.000000001
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$MR <- our$MR_edu
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((1+interval)/2, df = n - 1) * vec_sd / sqrt(n) # 1-interval or 1+interval is ==
  #error <- qnorm((1+interval)/2) * vec_sd / sqrt(n) # 1-interval or 1+interval is ==
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

#al = 1 - lev
#w = 1 - al/2
#1 - (1-lev)/2
#(  2 - 1 + lev )  /2
#https://math.stackexchange.com/questions/1350635/when-do-i-use-a-z-score-vs-a-t-score-for-confidence-intervals
#https://stats.stackexchange.com/questions/502482/confidence-intervals-95-we-do-z-qnorm-0-95-or-z-qnorm-0-975

LEVEL <- 0.99
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'MR'], LEVEL) })
colnames(ic_clones) <- unique(our$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

our$model_clone <- paste0(our$model, "_", our$clone)

n <- length(levels(as.factor(our$model_clone)))
# shape clones
our$time <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})
palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model_clone

save.image('mrca.Rdata')
ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ylab(yname)+
  geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), legend.position="none", axis.title.y=element_text(size=15))+scale_color_manual(values=pal)+scale_shape_manual(values=c(18,20))

ggsave(outfile)