lm_vafs  <- snakemake@input[['LM']]
pr_vafs  <- snakemake@input[['PR']]

outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))
library(ggplot2)
library(ggpubr)
load(theme)


load_data <- function(file) {
  data <- read.table(gzfile(file), sep="\t", header=TRUE)
  afcolumn <- colnames(data)[grepl('CRC1599', colnames(data))]
  af <- data[,afcolumn, drop=TRUE]
  return(data.frame(VAF=af))
}

dP <- load_data(pr_vafs)
dM <- load_data(lm_vafs)
x_breaks<-guess_ticks(c(0,1))
y_breaks<-guess_ticks(c(0,21))
print(dM$VAF)
print(x_breaks)
#hM <- hist(afM, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col="#ff9900", main="")
#hP <- hist(afP, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col="#ffff00", main="", ylim=c(0,25))

Mp <- ggplot(data=dM, aes(x=VAF))+geom_histogram( fill="#ff9900", color="black", size=0.2,binwidth=(max(dM$VAF)-min(dM$VAF))/49)+#ylim(0, 21)+ #bins=50, ,
  ylab('SNVs')+unmute_theme+ggtitle('CRC1599LMX')+
  geom_vline(xintercept=0.24, linetype="dashed", color = "darkgrey", size=0.2)+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
Pp <- ggplot(data=dP, aes(x=VAF))+geom_histogram(fill="#ffcc33", color="black", size=0.2,binwidth=(max(dP$VAF)-min(dP$VAF))/49)+#ylim(0, 21)+ ,boundary=max(dP$VAF)/50
  ylab('SNVs')+unmute_theme+ggtitle('CRC1599PRX')+
  geom_vline(xintercept=0.24, linetype="dashed", color = "darkgrey", size=0.2)+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))

p <- ggarrange(Pp, Mp,  common.legend = TRUE)
pdf('fig_5b_VAFhistogram.pdf')
print(p)
graphics.off()
ggsave(outplot, plot=p, width=120, height=60, units="mm")



save.image(paste0(outplot, '.Rdata'))