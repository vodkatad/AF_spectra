#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
infile1 <- args[1]
infile2 <- args[2]
outputname <- args[3]
colors <- args[4]
ggtheme <- args[5]
load(ggtheme)

outputplot <- paste0(outputname, '.pdf')
outputstats <- paste0(outputname, '.log')

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model_clone
load(ggtheme)

names <- strsplit(outputname, "_vs_")
name1 <- names[[1]][1]
name2 <- names[[1]][2]

data1 <- read.table(infile1, sep="\t", stringsAsFactors = FALSE)
data2 <- read.table(infile2, sep="\t", stringsAsFactors = FALSE)
colnames(data1) <- c('sample', name1)
colnames(data2) <- c('sample', name2)

m <- merge(data1, data2, by="sample")

m[, name1] <- m[, name1] / 0.000000001
m[, name2] <- m[, name2] / 0.000000001
m$model <- sapply(m$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
m$clone <- sapply(m$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
m$clone2 <- sapply(m$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
m$model_clone <- paste0(m$model, "_", m$clone)
m$time <- sapply(m$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

models <- c('CRC1307',  'CRC1502')
sink(outputstats)
print('Pearson on single clones estimates')
print(models)
sink()
for (mi in models) {
  sink(outputstats, append=TRUE)
  onem <- m[m$model == mi,]
  print(mi)
  print(cor.test(onem[, name1], onem[, name2]))
  sink()
}

#ggplot(data=m, aes(x=platy, y=noplaty, color=model_clone))+geom_point()+
#  unmute_theme+scale_color_manual(values=pal)
#ggplot(data=m, aes(x=cov1x, y=cov20x, color=model_clone))+geom_point()+
#  unmute_theme+scale_color_manual(values=pal)



prepare_plot <- function(our) {
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
  ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'MR'], LEVEL) })
  colnames(ic_clones) <- unique(our$model)
  pdata <- as.data.frame(t(ic_clones))
  pdata$model <- rownames(pdata)
  return(pdata)
}

our1 <- m
our1[,name2] <- NULL
colnames(our1)[colnames(our1) == name1] <- 'MR'

our2 <- m
our2[,name1] <- NULL
colnames(our2)[colnames(our2) == name2] <- 'MR'

pdata1 <- prepare_plot(our1)
pdata2 <- prepare_plot(our2)
pdata1$class <- name1
pdata2$class <- name2

our1$class <- name1
our2$class <- name2

our <- rbind(our1, our2)
pdata <- rbind(pdata1, pdata2)

ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
  geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=pal)+scale_shape_manual(values=c(18,20))+facet_wrap(~class)

ggsave(outputplot)

sink(outputstats, append=TRUE)
print(pdata1$mean / pdata2$mean)
sink()

save.image('pippo.RData')