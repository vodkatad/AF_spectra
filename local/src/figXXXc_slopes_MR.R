williams_f  <- snakemake@input[['williams']]
MR_f  <- snakemake@input[['MR']]
log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
WES_bp <- as.numeric(snakemake@params[['WESbp']]) # should correct for CN? mh
save.image(paste0(outplot, '.Rdata'))
library(ggplot2)
library(ggpubr)
library(reshape)
load(theme)

wdata <- read.table(williams_f, sep="\t", header=TRUE)
wdata$type <- substr(rownames(wdata), 0, 9)
MR <- read.table(MR_f, sep="\t", header=TRUE)

MR <- MR[grepl('CRC1599', MR$model),]
MR$type <- substr(MR$model, 0, 9)
wdata <- wdata[grepl("CRC1599", wdata$type),]
wdata <- wdata[wdata$r > 0.90 & wdata$subcl > 10,]

m <- merge(wdata, MR, by="type")
m$MR_b <- m$intercept / WES_bp

m$MR_MA <- m$mean 
m$MR_emp <- m$MR_b

m$type <- factor(m$model, levels= c('CRC1599PR', 'CRC1599LM'))

pdl <- melt(m, id='type', measure.vars = c('MR_MA', 'MR_emp'))

guess_neg_ticks <- function (values, nticks = 5, fixed_max = 0, fixed_min = NULL) 
{
  vmin <- min(values)
  if (is.null(fixed_min)) {
    round_min <- floor(vmin)
  }
  else {
    round_min <- fixed_min
  }
  my_breaks <- seq(round_min, fixed_max, length.out = nticks)
  return(my_breaks)
}

pdl$log <- log10(pdl$value)
y_breaks<- guess_neg_ticks(pdl$log)

p <- ggplot(data=pdl, aes(x=variable, y=log, fill=type)) +
  geom_col(position="dodge")+
  scale_fill_manual(values=c('#ffcc33', '#ff9900'))+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),0), expand = c(0, 0))+
     unmute_theme
  #   ylab('MR empiric')+xlab('MR MA')+unmute_theme+theme(legend.position="none")+
  #   scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  #   scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)), expand = c(0, 0))

# x_breaks<-guess_ticks(MR$mean)
# 
# p <- ggplot(data=m, aes(x=mean, y=MR_b, color=type)) +
#   geom_point(size=1)+
#   #geom_errorbar(data=pdata, aes(x=model, y=lower, ymin=lower, ymax=upper,group=NULL, color=NULL), width=.1, size=.2)+
#   scale_color_manual(values=c('#ffcc33', '#ff9900'))+
#   #scale_y_continuous(breaks = seq(0, 70, by = 10))+
#   ylab('MR empiric')+xlab('MR MA')+unmute_theme+theme(legend.position="none")+
#   scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
#   scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)), expand = c(0, 0))
  


ggsave(outplot, plot=p, width=60, height=60, units="mm")


pp <- p + theme(legend.position= "none")
ggsave(paste0('nolegend_', outplot), plot=pp, width=60, height=60, units="mm")

save.image(paste0(outplot, '.Rdata'))

q(0)
###
setwd('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/')
load('fig_XXXc_slopes_MR.svg.Rdata')

### tentativo 1
MR_all <- read.table('../datasetV2/MR_edu_SNV', sep="\t", stringsAsFactors = FALSE)
colnames(MR_all) <- c('sample', 'MR')
MR_all <- MR_all[grepl("CRC1599", MR_all$sample),]
MR_all$time <- sapply(MR_all$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})
MR_all$type <- substr(MR_all$sample, 0, 9)
MR_all <- MR_all[MR_all$time == 1,]

wdata$MR_b <- wdata$intercept / WES_bp

mm <- merge(MR_all, wdata, by="type")

mm$MR <- mm$MR / 0.000000001
mm$MR_b <- mm$MR_b / 0.000000001

y_breaks<- guess_ticks(mm$MR_b, fixed_max=635)
x_breaks<- guess_ticks(mm$MR)

mm$type <- factor(mm$type, levels= c('CRC1599PR', 'CRC1599LM'))

p <- ggplot(data=mm, aes(x=MR, y=MR_b, color=type)) +
  geom_point(size=1)+
  scale_color_manual(values=c('#ffcc33', '#ff9900'))+
  scale_y_continuous(breaks=y_breaks, limits=c(0, max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0, max(x_breaks)), expand = c(0, 0))+
  unmute_theme+ theme(legend.position= "none")

ggsave('mah.svg', plot=p, width=60, height=60, units="mm")

### tentativo 2
ttable <- m[, c('type', 'MR_MA', 'MR_emp')]
colnames(ttable) <- c('model', "MA_expr", 'MA_is')
ttable$foldchange <- ttable$MA_is /  ttable$MA_expr
ttable$model <- NULL
ttable$foldchange <- round(ttable$foldchange, digits = 2)
ttable <- rbind(ttable, c(ttable[1,1]/ttable[1,2]))
write.table(ttable, file="mah.tsv", sep="\t", quote=FALSE, row.names = FALSE)
