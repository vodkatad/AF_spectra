library(ggplot2)

ratios_f  <- snakemake@input[['sign']]
log_f <- snakemake@log[['log']]
outplot1 <- snakemake@output[['plot1']]

theme <- snakemake@input[['theme']]
load(theme)
save.image(paste0(outplot1, '.Rdata'))

d <- read.table(gzfile(ratios_f), sep="\t", header=T)

d$lr_cl <- log((d$Signature.8_cl/d$Signature.1_cl)+1)
d$lr_sl <- log((d$Signature.8_subcl/d$Signature.1_subcl)+1) 

d_nona <- d[!is.infinite(d$lr_cl) & !is.infinite(d$lr_sl) & !is.na(d$lr_cl) & !is.na(d$lr_sl),]

pd <- data.frame(ratio=c(d_nona$lr_cl, d_nona$lr_sl), class=c(rep('Clonal', nrow(d_nona)), rep('Subclonal', nrow(d_nona))))


sink(log_f)
print('Starting pr:')
nrow(d)
print('Removed NA:')
nrow(d_nona)
sink()

save.image(paste0(outplot1, '.Rdata'))

y_breaks <- guess_ticks(pd$ratio)

p <- ggplot(data=pd, aes(y=ratio, x=class))+geom_boxplot(outlier.size=1)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  unmute_theme+ylab('log(ratio+1)')+xlab('SNVs')


w <- wilcox.test(d$ratio_subcl, d$ratio_cl, alternative="greater", paired=TRUE)
tt <- wilcox.test(d$lr_sl, d$lr_cl, alternative="greater", paired=TRUE)

sink(log_f, append=T)
w
w$p.value
nrow(d)
tt
tt$p.value
sink()

save.image(paste0(outplot1, '.Rdata'))

ggsave(outplot1, plot=p, width=89, height=89, units="mm")


q('no') # slides

textSize <- 10
largerSize <- textSize + 2

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
unmute_theme <- theme(
  text = element_text(size = textSize, family='sans'),
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(color = "black"),
  panel.background = element_blank()
)
sl <- ggplot(data=pd, aes(y=ratio, x=class))+geom_boxplot(outlier.size=1)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  unmute_theme+ylab('log(ratio+1)')+xlab('SNVs')



ggsave(sl, file="~/slide1.pdf", width=89, height=89, units="mm")
