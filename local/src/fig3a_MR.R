MR_f  <- snakemake@input[['MR']]
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

colnames(our) <- c('sample','MR_edu')

our$MR <- our$MR_edu / 0.000000001
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)
our$time <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

# we keep only models with at least one t2
models_2nd_round <- unique(our[our$time==2, 'model'])
our <- our[our$model %in% models_2nd_round,]

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

# we compute avg on model_time
our$model_time <- paste0(our$model, "_", our$time)
LEVEL <- 0.99
ic_clones <- sapply(unique(our$model_time), function(x) { confidence_interval(our[our$model_time==x,'MR'], LEVEL) })
colnames(ic_clones) <- unique(our$model_time)
pdata <- as.data.frame(t(ic_clones))
pdata$model_time <- rownames(pdata)
pdata$model <- sapply(pdata$model_time, function(x) {y<-strsplit(x, '_')[[1]][1]; return(y[1])})
pdata$time <- sapply(pdata$model_time, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})


orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=TRUE)
orderdf$mean <- NULL
orderdf$lower <- NULL
orderdf$upper <- NULL
orderdf$model <- paste0(orderdf$model, ifelse(!grepl('\\d$', orderdf$model), '', ifelse(orderdf$model=="CRC0282", 'PR', 'LM')))
orderdf$x_ord <- seq(1, nrow(orderdf))
our <- merge(our, orderdf, by="model")
pdata <- merge(pdata, orderdf, by="model")
our$ord_x_time <- our$x_ord + 1/3*as.numeric(our$time)
pdata$ord_x_time <- pdata$x_ord + 1/3*as.numeric(pdata$time)


y_breaks <- guess_ticks(our$MR)
print(y_breaks)
our$model_time <- as.factor(our$model_time)
p <- ggplot(pdata, aes(x=reorder(model_time, ord_x_time), y=mean)) +  geom_point(stat="identity", shape=1, size=2) +
  geom_segment(aes(y=lower, yend=upper, x=reorder(model_time, ord_x_time), xend=reorder(model_time, ord_x_time)), size=0.6)+ylab('MR')+xlab('')+
  geom_point(data=our, aes(x=reorder(model_time, ord_x_time), y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.5))+
  scale_color_manual(values=pal)+unmute_theme+
  scale_y_continuous(breaks=y_breaks) +#, expand = c(0, 0))+
  theme(legend.position="right", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))

ggsave(outplot, plot=p, width=89, height=89, units="mm")
write.table(pdata, file=data_f, sep="\t", quote=FALSE)

save.image(paste0(outplot, '.Rdata'))
