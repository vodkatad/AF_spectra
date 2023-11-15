MR_f  <- snakemake@input[['MR']]
order_f  <- snakemake@input[['order']]
colors <- snakemake@input[['palette']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
data_f <- snakemake@output[['avgdata']]

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

our <- read.table(MR_f, sep="\t", header=FALSE, stringsAsFactors=FALSE)

colnames(our) <- c('sample','MR_edu')
our <- our[grepl('-1-', our$sample),] # only 1st round

our$MR <- our$MR_edu / 0.000000001 #/1000
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)

our$model <- paste0(our$model, ifelse(!grepl('\\d$', our$model), '', ifelse(our$model=="CRC0282", 'PR', 'LM')))

our <- our[our$model != "CRC0282PR",]

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
orderdf <- orderdf[orderdf$model != "CRC0282PR",]

our$model <- factor(our$model, levels=orderdf$model)
pdata$model <- factor(pdata$model, levels=orderdf$model)

y_breaks <- guess_ticks(our$MR, fixed_max=0.5)


#p <- ggplot() + 
#  geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.5))+
#  geom_point(data=pdata, aes(x=model, y=mean), stat="identity", shape=1, size=2) +
#  geom_segment(data=pdata, aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+ylab('MR')+xlab('')+
#  scale_color_manual(values=pal)+theme(legend.position="none", axis.text.x = element_blank(), axis.ticks = element_blank())+
#  scale_y_continuous(breaks=y_breaks, limits=min(y_breaks),max(y_breaks))+
#  unmute_theme
pdata$xmodel <- as.numeric(pdata$model)
#p <- ggplot() + 
  #geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.7))+
  #geom_segment(data=pdata, aes(x=xmodel-0.4, y=mean, yend=mean, xend=xmodel+0.4), size=0.3) +
  #geom_errorbar(data=pdata, aes(ymin=lower, ymax=upper, x=model), size=0.3, width=0.3)+ylab('MR (SNV/(Gbp*division))')+xlab('')+
  #scale_color_manual(values=pal)+
  #scale_y_continuous(breaks=y_breaks, expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  #unmute_theme+theme(legend.position="right", axis.text.x = element_blank(), 
                     #axis.ticks.x = element_blank(),
                     #legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))
#ratio_to_caperrorbars <- y_breaks[2]
#print('ratio')
#print(ratio_to_caperrorbars)
#pdf('fig_1b_MR.pdf')

p <- ggplot() + 
  geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.7))+
  geom_segment(data=pdata, aes(x=xmodel-0.2, yend=mean,y=mean,  xend=xmodel+0.2),size=.3) +
  geom_errorbar(data=pdata, aes(ymin=lower, ymax=upper, x=model), size=0.3, width=0.3)+ylab('MR [SNV/(Gbp*division)]')+xlab('PDTs')+
  scale_color_manual(values=pal)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   
#print(p)
#graphics.off()
#print(our)
ggsave(outplot, plot=p, width=89, height=89, units="mm")
write.table(pdata, file=data_f, sep="\t", quote=FALSE)
sink(log_f)
print('n clones')
print(nrow(our))

save.image(paste0(outplot, '.Rdata'))
