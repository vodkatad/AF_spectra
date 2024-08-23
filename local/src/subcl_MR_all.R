library(ggplot2)
library(ggpubr)

subclonal_f  <- snakemake@input[['subclonal']]
MR_f  <- snakemake@input[['MR']]
outplot <- snakemake@output[['plot']]
log_f <- snakemake@log[['log']]

colors <- snakemake@input[['palette']]
theme <- snakemake@input[['theme']]

load(theme)

save.image(paste0(outplot, '.Rdata'))

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model_clone

subclonal <- read.table(subclonal_f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(subclonal) <- c('sample', 'n')
MR <- read.table(MR_f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(MR) <- c('sample', 'MR')
MR$MR_edu <- MR$MR / 0.000000001
subclonal$n <- subclonal$n / 1000
merged <- merge(subclonal, MR, by='sample')
merged$smodel <- substr(merged$sample, 0, 7)

merged$model <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
merged$PDT <- paste0(merged$model, ifelse(!grepl('\\d$', merged$model), '', ifelse(merged$model=="CRC0282", 'PR', 'LM')))

merged$clone <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
merged$clone2 <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
merged$model_clone <- paste0(merged$model, "_", merged$clone)
merged$time <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

merged <- merged[merged$time != '2',]
sink(log_f)
print('N clones')
print(nrow(merged))
print('N models')
print(length(unique(merged$model)))
sink()

# average errors ############
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  #error <- qt((1+interval)/2, df = n - 1) * vec_sd / sqrt(n)
  error <- vec_sd
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

y_breaks <- guess_ticks(merged$n)
x_breaks <- guess_ticks(merged$MR_edu)
print(y_breaks)
print(x_breaks)


# textSize <- 15
# largerSize <- 20
# slidetheme <- theme(
#   text = element_text(size = textSize, family='sans'),
#   axis.title = element_text(size = largerSize),
#   axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
#   axis.text.y = element_text(size = textSize, color="black"),
#   plot.title = element_text(size = largerSize, hjust = 0.5),
#   legend.title = element_text(size=largerSize, hjust = 0.5),
#   legend.text = element_text(size=textSize),
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   axis.line = element_line(colour = "black"),
#   axis.ticks = element_line(color = "black"),
#   panel.background = element_blank()
# )


ratio_to_caperrorbars <- y_breaks[2] / x_breaks[2]
p <- ggplot() +
  geom_point(data=merged, aes(x=MR_edu, y=n, color=model_clone, fill=model_clone), stat="identity", size=0.4)+#, position=position_dodge(0.7))+
  unmute_theme+scale_color_manual(values=pal, guide="none")+#scale_shape_manual(values=c(18,23,20,19))+
  geom_point(data=pdata, aes(x=mean.x, y=mean.y), shape=1, size=1)+
  geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, ymin=lower.y, ymax=upper.y), width=.1, size=.2)+
  geom_errorbarh(data=pdata, aes(y=mean.y, xmin=lower.x, xmax=upper.x), height=.1*ratio_to_caperrorbars, size=.2)+
  scale_fill_manual(values=pal, guide="none")+
  xlab('MR [SNVs/(Gbp*division)]')+ylab('# subclonal SNVs/kbp')+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)), expand = c(0, 0))+
  guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))
pp <- p + theme(legend.position= "none")


# pearson pvalue
sink(log_f, append=TRUE)
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
  p <- cor.test(pdata$mean.x, pdata$mean.y, method=me)
  print(p)
  if (me == 'pearson') {
    print(p$estimate**2)
  }
  print(nrow(pdata))
}
sink()

ggsave(outplot, plot=pp, width=60, height=60, units="mm")
save.image(paste0(outplot, '.Rdata'))

