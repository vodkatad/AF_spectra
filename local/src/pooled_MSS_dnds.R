library(ggplot2)
t1 <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/dnds_vitroMSS_t1.tsv', sep="\t")
t2 <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/dnds_vitroMSS_t1t2.tsv', sep="\t")

d <- data.frame(x=c('T1', 'T1_T2'), estimate=c(t1['wall', 'mle'],t2['wall', 'mle']), upper=c(t1['wall', 'cihigh'],t2['wall', 'cihigh']), lower=c(t1['wall', 'cilow'],t2['wall', 'cilow']))

ggplot(d, aes(x=x, y=estimate)) +  geom_point(stat="identity", size=1.5) +
  geom_hline(yintercept=1,linetype=2,size=0.2)+
  geom_errorbar(aes(ymin=lower, ymax=upper, x=x, width=0.3), size=0.5, color='black')+ylab('dN/dS estimate')+xlab('MA round')+theme_bw()

### and gecip
theme <- '/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_5.Rdata'

library(ggplot2)
load(theme)


d <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/ratios_8_1.tsv', sep="\t", header=T)
dnona <- d[!is.na(d$ratio),]
dnona$x <- ifelse(dnona$type =="PRIMARY", 'PRs', 'LMs')
dnona$x <- factor(dnona$x, levels=c('PRs', 'LMs'))

y_breaks <- guess_ticks(dnona$ratio)
p <- ggplot(data=dnona, aes(x=x, y=ratio))+geom_violin()+geom_boxplot(width=0.1,outlier.size = 0.5)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  unmute_theme+xlab('Type')+ylab('SBS8/SBS1')
#+theme(legend.position="none", axis.text.x = element_blank(), 
#                     axis.ticks.x = element_blank(),
#                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))


ggsave('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/manual_1.svg', plot=p, width=89, height=89, units="mm")


nsample <- nrow(dnona[dnona$type=="METASTASES", ])

pri <- dnona[dnona$type!="METASTASES", ]
met <- dnona[dnona$type=="METASTASES", ]

onesample <- function(mpri, mmet, n) {
  smpri <-mpri[sample(rownames(mpri), n),]
  #wi <- wilcox.test(mmet[, 'ratio'], smpri[, 'ratio'], side="greater")
  #return(wi$p.value)
  return(median(mmet$ratio) - median(smpri$ratio))
}

set.seed(42)
pivs <- replicate(10000, onesample(pri, met, nsample))

#https://stackoverflow.com/questions/77589770/using-limits-in-scale-x-continuous-with-geom-histogram-removes-values-even-when
#hist(-log10(pivs), breaks=50); abline(v=-log10(0.05), col="red")
hist(pivs, breaks=50); abline(v=0, col="red")
y_breaks <- guess_ticks(c(0, 1000))

pd <- data.frame(n=pivs)
p<-ggplot(data=pd, aes(x=n))+geom_histogram(bins=40, fill="white", color="black")+geom_vline(xintercept=0, color='red')+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+coord_cartesian(xlim = c(-0.4, 0.4))+
  unmute_theme+ylab('# Simulations')+xlab('Median LMs - Median PRs')
ggsave('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/manual_2.svg', plot=p, width=89, height=89, units="mm")


table(pivs>0.05)
63/10000*100
