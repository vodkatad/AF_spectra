edu_f  <- snakemake@input[['edu']]
order_f  <- snakemake@input[['order']]
log_f  <- snakemake@log[['log']]
colors <- snakemake@input[['palette']]

outplot <- snakemake@output[['plot']]

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model


our <- read.table(edu_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)


colnames(our) <- c('Campione', 'sample', 'conc', 'ncell', 'cell_date', 'analysis_date', 'EDU', 'n_div', 'T')
# right now keep only T1 
our <- our[our$T == 1,]
#our$well <- ifelse(grepl('100', our$Campione), '96_wells', '12_wells')
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})

our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)

our$model <- paste0(our$model, ifelse(!grepl('\\d$', our$model), '', ifelse(our$model=="CRC0282", 'PR', 'LM')))

#our$model_well <- paste0(our$model, "_", our$well)
#our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
#our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})

sink(log_f)
print(nrow(our))
print(length(unique(our$model_clone)))
print(table(our$model))
sink()

#figura1
#colnames(our) <- c('sample','MR_edu')
#our <- our[grepl('-1-', our$sample),] # only 1st round

#our$MR <- our$MR_edu/ 0.000000001 #/1000


#fine figura 1

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
  if (is.na(error)) {
    result <- c("lower" = vec_mean, "upper" = vec_mean, "mean" = vec_mean)
  }
  return(result)
}

# we compute avg on model_time
LEVEL <- 0.99
#ic_clones <- sapply(unique(our$model_well), function(x) { confidence_interval(our[our$model_well==x,'EDU'], LEVEL) })
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'EDU'], LEVEL) })
colnames(ic_clones) <- unique(our$model)#colnames(ic_clones) <- unique(our$model_well)
pdata <- as.data.frame(t(ic_clones))
#pdata$model_well <- rownames(pdata)
pdata$model <- rownames(pdata) #pdata$model <- sapply(pdata$model_well, function(x) {y<-strsplit(x, '_')[[1]][1]; return(y[1])})
#pdata$well <- sapply(pdata$model_well, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})

orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=TRUE)
orderdf$model <- paste0(orderdf$model, ifelse(!grepl('\\d$', orderdf$model), '', ifelse(orderdf$model=="CRC0282", 'PR', 'LM')))
our$model <- factor(our$model, levels=orderdf$model)
pdata$model <- factor(pdata$model, levels=orderdf$model)
pdata$xmodel <- as.numeric(pdata$model)
#orderdf$mean <- NULL
#orderdf$lower <- NULL
#orderdf$upper <- NULL
#orderdf$model <- paste0(orderdf$model, ifelse(!grepl('\\d$', orderdf$model), '', ifelse(orderdf$model=="CRC0282", 'PR', 'LM')))
#orderdf <- orderdf[orderdf$model %in% unique(pdata$model),, drop=FALSE]
#orderdf$x_ord <- seq(1, nrow(orderdf))
#our <- merge(our, orderdf, by="model")
#pdata <- merge(pdata, orderdf, by="model")
#our$ord_x_time <- our$x_ord + 0.4*ifelse(our$well == '96_wells', 1, 2)
#pdata$ord_x_time <- pdata$x_ord +0.4*ifelse(pdata$well == 96, 1, 2)


y_breaks <- guess_ticks(our$EDU,fixed_max=50)

##IANG
#names(pal) <- gsub('-','_', names(pal))
#our$model_well <- as.factor(our$model_well)
p <- ggplot(pdata, aes(x=model, y=mean))+  #geom_point(stat="identity", shape=1, size=2) +
  geom_point(data=our, aes(x= model, y=EDU, color=model_clone), stat="identity", size=2, shape=18,position=position_dodge(0.3))+#, 
  geom_segment(data=pdata, aes(x= xmodel-0.2, yend=mean,y=mean,  xend= xmodel+0.2),size=.3)+
  geom_errorbar(aes(ymin=lower, ymax=upper, x= xmodel), size=0.3, width=0.3)+ylab('% cells EDU+')+xlab('PDTs')+#xmin=reorder(model_time, ord_x_time), xmax=reorder(model_time, ord_x_time)
  scale_color_manual(values=pal)+unmute_theme+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0)) +#, expand = c(0, 0))+
  theme(legend.position="none", 
                    axis.ticks.x = element_blank(),axis.text.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))
## why -0.05 ?? TODO FIXME
#axis.text.x = element_blank(),
ggsave(outplot, plot=p, width=89, height=89, units="mm")

save.image(paste0(outplot, '.Rdata'))