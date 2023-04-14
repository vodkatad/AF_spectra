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

#hM <- hist(afM, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col="#ff9900", main="")
#hP <- hist(afP, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col="#ffff00", main="", ylim=c(0,25))

Mp <- ggplot(data=dM, aes(x=VAF))+geom_histogram(bins=50, fill="#ff9900", color="black", size=0.2)+ylim(0, 21)+
  ylab('')+unmute_theme+ggtitle('CRC1599LMX')+
  geom_vline(xintercept=0.24, linetype="dashed", color = "darkgrey", size=0.2)
Pp <- ggplot(data=dP, aes(x=VAF))+geom_histogram(bins=50, fill="#ffff00", color="black", size=0.2)+ylim(0, 21)+
  ylab('# of alterations')+unmute_theme+ggtitle('CRC1599PRX')+
  geom_vline(xintercept=0.24, linetype="dashed", color = "darkgrey", size=0.2)

p <- ggarrange(Pp, Mp,  common.legend = TRUE)

ggsave(outplot, plot=p, width=89, height=56, units="mm")



save.image(paste0(outplot, '.Rdata'))