library(ggplot2)

ratios_f  <- snakemake@input[['ratios']]
log_f <- snakemake@log[['log']]
outplot1 <- snakemake@output[['plot1']]
nmonte <- as.numeric(snakemake@params[['random']])
nsim <- as.numeric(snakemake@params[['nsim']])

theme <- snakemake@input[['theme']]
load(theme)
save.image(paste0(outplot1, '.Rdata'))

d <- read.table(ratios_f, sep="\t", header=T)
dnona <- d[!is.na(d$ratio),]
dnona$x <- ifelse(dnona$type =="PRIMARY", 'PRs', 'LMs')
dnona$x <- factor(dnona$x, levels=c('PRs', 'LMs'))

nsample <- nrow(dnona[dnona$type=="METASTASES", ]) # avoid picking always all mets
sink(log_f)
print('Nrand=')
print(nmonte)
print(nsample)
sink()
nsample <- nmonte

pri <- dnona[dnona$type!="METASTASES", ]
met <- dnona[dnona$type=="METASTASES", ]

sink(log_f, append=TRUE)
print('Npri=')
print(nrow(pri))
print('Nmet=')
print(nrow(met))
print('NisNA=')
print(table(is.na(d$ratio)))
print(table(d$type))
print(table(dnona$type))
sink()
onesample <- function(mpri, mmet, n) {
  smpri <- mpri[sample(rownames(mpri), n),]
  smmet <- mmet[sample(rownames(mmet), n),]
  #wi <- wilcox.test(mmet[, 'ratio'], smpri[, 'ratio'], side="greater")
  #return(wi$p.value)
  return(median(smmet$ratio) - median(smpri$ratio))
}

set.seed(42)
pivs <- replicate(nsim, onesample(pri, met, nsample))

#https://stackoverflow.com/questions/77589770/using-limits-in-scale-x-continuous-with-geom-histogram-removes-values-even-when
#hist(-log10(pivs), breaks=50); abline(v=-log10(0.05), col="red")
#hist(pivs, breaks=50); abline(v=0, col="red")
y_breaks <- guess_ticks(c(0, 1000))

pd <- data.frame(n=pivs)
p <- ggplot(data=pd, aes(x=n))+geom_histogram(bins=40, fill="white", color="black")+geom_vline(xintercept=0, color='red')+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+coord_cartesian(xlim = c(-0.4, 0.4))+
  unmute_theme+ylab('# Simulations')+xlab('Median LMs - Median PRs')

ggsave(outplot1, plot=p, width=89, height=89, units="mm")


dd <- as.data.frame(table(pivs>0))
sink(log_f, append=TRUE)
dd
dd[dd$Var1==T, 'Freq']/nsim
dd[dd$Var1==F, 'Freq']/nsim
sink()
save.image(paste0(outplot1, '.Rdata'))

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
sl <- ggplot(data=pd, aes(x=n))+geom_histogram(bins=40, fill="white", color="black")+geom_vline(xintercept=0, color='red')+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+coord_cartesian(xlim = c(-0.4, 0.4))+
  unmute_theme+ylab('# Simulations')+xlab('Median LMs - Median PRs')


ggsave(sl, file="~/slide1.pdf", width=89, height=89, units="mm")
