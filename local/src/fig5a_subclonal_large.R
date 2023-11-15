subclonal_f  <- snakemake@input[['subclonal']]
MR_f  <- snakemake@input[['MR']]
colors <- snakemake@input[['colors']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
outtsv <- snakemake@output[['tsv']]
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
subclonal$n <- subclonal$n / 1000
merged <- merge(subclonal, MR, by='sample')
merged$smodel <- substr(merged$sample, 0, 7)
wanted <- c('CRC1599', 'CRC1307', 'CRC1078')
merged <- merged[merged$smodel %in% wanted,]
sink(log_f)
print('N clones')
print(nrow(merged))
sink()

merged$model <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
merged$PDT <- paste0(merged$model, ifelse(grepl('\\d$', merged$model), 'LM', ''))
merged$PDT <- factor(merged$PDT, levels=c('CRC1078LM','CRC1307LM','CRC1599PR','CRC1599LM'))

merged$clone <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
merged$clone2 <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
merged$model_clone <- paste0(merged$model, "_", merged$clone)
merged$time <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

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

sink(log_f, append=TRUE)
# fold change + wilcox CRC1307 CRC1078
m1 <- "CRC1307"
m2 <- "CRC1078"
print(paste0(m1, ' vs ', m2))
op <- pdataMR[order(pdataMR$mean),]
print(as.character(op[nrow(op), 'model']))
print(as.character(op[nrow(op)-1, 'model']))
print(pdataMR[pdataMR$model==m1, 'mean']/pdataMR[pdataMR$model==m2, 'mean'])
print(wilcox.test(merged[merged$model==m1, 'MR'], merged[merged$model==m2, 'MR']))

# fold change + wilcox CRC1599PR CRC1599LM
m1 <- "CRC1599LM"
m2 <- "CRC1599PR"
print(paste0(m1, ' vs ', m2))
print(as.character(op[1, 'model']))
print(as.character(op[nrow(op)-1, 'model']))
print(pdataMR[pdataMR$model==m1, 'mean']/pdataMR[pdataMR$model==m2, 'mean'])
print(wilcox.test(merged[merged$model==m1, 'MR'], merged[merged$model==m2, 'MR']))
sink()

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

ratio_to_caperrorbars <- y_breaks[2] / x_breaks[2]
p <- ggplot() +
   geom_point(data=merged, aes(x=MR_edu, y=n, color=model_clone, shape=PDT, fill=model_clone), stat="identity", size=1, position=position_dodge(0.7))+
   unmute_theme+scale_color_manual(values=pal, guide="none")+scale_shape_manual(values=c(24,25,22,23))+#scale_shape_manual(values=c(18,23,20,19))+
   geom_point(data=pdata, aes(x=mean.x, y=mean.y), shape=1, size=2)+
   geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, ymin=lower.y, ymax=upper.y), width=.1, size=.2)+
   geom_errorbarh(data=pdata, aes(y=mean.y, xmin=lower.x, xmax=upper.x), height=.1*ratio_to_caperrorbars, size=.2)+
   scale_fill_manual(values=pal, guide="none")+
   xlab('MR [SNVs/(Gbp*division)]')+ylab('# subclonal SNVs')+
   scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
   scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)), expand = c(0, 0))+
   theme(legend.position="none",  legend.spacing.y = unit(0.15, "mm")) + 
   guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))


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
  print(cor.test(pdata$mean.x, pdata$mean.y, method=me))
  print(nrow(pdata))
}
sink()
pdf('fig_5a_subclonal.pdf')
print(p)
graphics.off()
ggsave(outplot, plot=p, width=89, height=89, units="mm")

pp <- p + theme(legend.position= "none")
ggsave(paste0('nolegend_', outplot), plot=pp, width=60, height=60, units="mm")

res <- merged[, c('sample', 'n')]
colnames(res) <- c('clone_id', 'n_subclonal')
write.table(res, file=outtsv, sep="\t", quote=FALSE, row.names=FALSE)

save.image(paste0(outplot, '.Rdata'))