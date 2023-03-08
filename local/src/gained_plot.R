#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
colors <- args[3]
vivi <- args[4]
norm <- args[5]

our <- read.table(infile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(our) <- c('sample','gained')
ylab <- 'N'
if (norm == "norm") {
  our$gained <- our$gained * 1000000
  ylab <- 'mut/Mpb'
}
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$vivi <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})


confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

LEVEL <- 0.99
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'gained'], LEVEL) })
colnames(ic_clones) <- unique(our$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

if (!is.na(vivi) && vivi == "vivo") {
  our$model_clone <- paste0(our$model, "_", our$clone, "_", our$vivi)
  scale_shape <- c(18,20,17)
  legend_position <- 'right'
} else {
  our$model_clone <- paste0(our$model, "_", our$clone)
  scale_shape <- c(18,20)
  legend_position <- 'none'
}

n <- length(levels(as.factor(our$model_clone)))

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model_clone


ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)




our$time <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})
if (n <= length(pal)) {
ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('Gained muts')+ylab(ylab)+xlab('')+
  geom_point(data=our, aes(x=model, y=gained, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=pal)+scale_shape_manual(values=scale_shape)
} else {
ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('Gained muts')+ylab(ylab)+xlab('')+
  geom_point(data=our, aes(x=model, y=gained, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_shape_manual(values=scale_shape)
}
ggsave(outfile)
save.image(paste0(outfile, '.Rdata'))