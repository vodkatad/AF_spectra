library(ggplot2)
library(reshape)
d <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/MA_treats_Vale/CTG_CRC1430_CRC1620_Afa.tsv', sep="\t", header=T)
colnames(d) <- gsub('AFA.', '', colnames(d), fixed=T)
colnames(d) <- gsub('.nM', '', colnames(d), fixed=T)
colnames(d) <- gsub('nM', '', colnames(d), fixed=T)

ave <- d[d$measure == 'ave',]
devstd <- d[d$measure == 'devstd',]
ave$measure <- NULL
devstd$measure <- NULL

lave <- melt(ave)
ldev <- melt(devstd)

lave$variable <- as.numeric(as.character(lave$variable))
ldev$variable <- as.numeric(as.character(ldev$variable))
colnames(lave) <- c('model', 'afa', 'ave')
colnames(ldev) <- c('dmodel', 'dafa', 'stdev')
pd <- cbind(lave, ldev)
ggplot(data=pd, aes(x=afa,y=ave))+geom_point()+theme_bw(base_size=15)+facet_wrap(~model)+
ylab('Relative cell number (T/NT)')+xlab('Afatinib (nM)')+scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1))+
geom_errorbar(aes(ymin = ave-stdev, ymax = ave+stdev),width = 0.01)+scale_x_log10(breaks=c(3, 6, 12.5, 25, 200))+
  geom_smooth(method="loess", se =FALSE)

ggplot(data=pd, aes(x=afa,y=ave))+geom_point()+theme_bw(base_size=15)+facet_wrap(~model)+
  ylab('Relative cell number (T/NT)')+xlab('Afatinib (nM)')+scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1))+
  geom_errorbar(aes(ymin = ave-stdev, ymax = ave+stdev),width = 0.01)+scale_x_log10(breaks=c(3, 6, 12.5, 25, 200))
  