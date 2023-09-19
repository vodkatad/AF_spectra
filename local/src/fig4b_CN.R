MR_f  <- snakemake@input[['MR']]
order_f  <- snakemake@input[['order']]
colors <- snakemake@input[['palette']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
data_f <- snakemake@output[['avgdata']]

rule <- snakemake@rule

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
library(ggpubr)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

our <- read.table(MR_f, sep="\t", header=FALSE, stringsAsFactors=FALSE)

colnames(our) <- c('sample','CN')
#our <- our[grepl('-1-', our$sample),] # only 1st round

our$MR <- our$CN
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)

our$model <- paste0(our$model, ifelse(!grepl('\\d$', our$model), '', ifelse(our$model=="CRC0282", 'PR', 'LM')))

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- vec_sd
  
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

LEVEL <- 0.99
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'MR'], LEVEL) })
colnames(ic_clones) <- unique(our$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=TRUE)
orderdf$model <- paste0(orderdf$model, ifelse(!grepl('\\d$', orderdf$model), '', ifelse(orderdf$model=="CRC0282", 'PR', 'LM')))
our$model <- factor(our$model, levels=orderdf$model)
pdata$model <- factor(pdata$model, levels=orderdf$model)
if (rule == "figure_4b_mr") {
  y_breaks <- guess_ticks(our$MR+1,fixed_max=60)
  ytitle <- 'MEDICC2 events'
} else if (rule == "figure_4_support_MR_also2nd") {
  y_breaks <- guess_ticks(our$MR, fixed_max=max(our$MR))
  ytitle <- 'MR'
} else {
  y_breaks <- guess_ticks(our$MR)
  ytitle <- 'tree distances (in bis n gained)'
}
pdata$xmodel <- as.numeric(pdata$model)

pdf('fig_4b_CN.pdf')
p <-  ggplot()+
#ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=2) +
  geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.5))+
  geom_errorbar(data=pdata,aes(ymin=lower, ymax=upper, x=model), size=0.3, width=0.3)+ylab(ytitle)+xlab('PDTs')+#, xmin=model
  geom_segment(data=pdata, aes(x=xmodel-0.2, yend=mean,y=mean,  xend=xmodel+0.2),size=.3) +
  scale_color_manual(values=pal)+unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm"))+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))
print(p)
graphics.off()
ggsave(outplot, plot=p, width=89, height=89, units="mm")
write.table(pdata, file=data_f, sep="\t", quote=FALSE)

sink(log_f)
print('n clones')
print(nrow(our))
sink()
save.image(paste0(outplot, '.Rdata'))
