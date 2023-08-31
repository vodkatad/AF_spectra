MR_f  <- snakemake@input[['N']]
order_f  <- snakemake@input[['order']]
colors <- snakemake@input[['palette']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
data_f <- snakemake@output[['avgdata']]

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
library(ggpubr)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

our <- read.table(MR_f, sep="\t", header=FALSE, stringsAsFactors=FALSE)

colnames(our) <- c('sample','MR')

our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)
our$time <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

# we keep only models with at least one t2
models_2nd_round <- unique(our[our$time==2, 'model'])
our <- our[our$model %in% models_2nd_round,]

our$model <- paste0(our$model, ifelse(!grepl('\\d$', our$model), '', ifelse(our$model=="CRC0282", 'PR', 'LM')))
# we compute avg on model_time
our$model_time <- paste0(our$model, "_", our$time)

orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=TRUE)
orderdf$mean <- NULL
orderdf$lower <- NULL
orderdf$upper <- NULL
orderdf$model <- paste0(orderdf$model, ifelse(!grepl('\\d$', orderdf$model), '', ifelse(orderdf$model=="CRC0282", 'PR', 'LM')))
orderdf <- orderdf[orderdf$model %in% unique(our$model),, drop=FALSE]
orderdf$x_ord <- seq(1, nrow(orderdf))
our <- merge(our, orderdf, by="model")
our$ord_x_time <- our$x_ord +as.numeric(our$time)
print(our$ord_x_time)
print('===========================0')
print(our$x_ord)
print('===========================0')
print(our$time)
print('================================')
print(orderdf)

our <- our[our$model=="CRC1502LM",]
print(our$model_clone)
y_breaks <- guess_ticks(our$MR, fixed_max=4800, fixed_min=2000)
our$model_time <- as.factor(our$model_time)
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
ic_clones <- sapply(unique(our$ord_x_time), function(x) {confidence_interval(our[our$ord_x_time==x,'MR'], LEVEL) })
colnames(ic_clones) <- unique(our$ord_x_time)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)
pdata$model<-as.numeric(pdata$model)
print(pdata)

p <- ggplot()+  #geom_point(stat="identity", shape=1, size=2) 
  #geom_hline(yintercept=1,linetype=2,size=0.2)+
  geom_point(data=our, aes(x= ord_x_time, y=MR, color=model_clone), stat="identity", size=2, shape=18,position=position_dodge(0.1))+#, 
  geom_segment(data=pdata, aes(x=model-0.2, yend=mean,y=mean,  xend=model+0.2),size=.3) +
  geom_errorbar(data=pdata, aes(ymin=lower, ymax=upper, x=model), size=0.3, width=0.3)+
  ylab('# SNVs')+xlab('PDTs')+#xmin=reorder(model_time, ord_x_time), xmax=reorder(model_time, ord_x_time)
  scale_color_manual(values=pal)+unmute_theme+
  scale_y_continuous(breaks=y_breaks,limits=c(2000, max(y_breaks)),expand = c(0, 0)) +#, expand = c(0, 0))+
  theme(legend.position="right", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))
#axis.text.x = element_blank(),
ggsave(outplot, plot=p, width=50, height=89, units="mm")

save.image(paste0(outplot, '.Rdata'))
 #ggplot(pdata, aes(x=reorder(model_time, ord_x_time), y=mean)) +  #geom_point(stat="identity", shape=1, size=2) +
  #geom_segment(data=pdata, aes(x=reorder(model_time, ord_x_time), yend=mean,y=mean,  xend=reorder(model_time, ord_x_time))),size=.2)+
  #geom_point(data=our, aes(x=reorder(model_time, ord_x_time), y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.5))+
  #geom_errorbar(aes(ymin=lower, ymax=upper, x=as.numeric(reorder(model_time, ord_x_time))), size=0.5)+ylab('MR')+xlab('')+#xmin=reorder(model_time, ord_x_time), xmax=reorder(model_time, ord_x_time)
  #scale_color_manual(values=pal)+unmute_theme+
  #scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0)) +#, expand = c(0, 0))+
  #theme(legend.position="right", axis.text.x = element_blank(), 
                     #axis.ticks.x = element_blank(),
                     #legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))