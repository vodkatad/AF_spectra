library(ggplot2)

ratios_f  <- snakemake@input[['ratios']]
log_f <- snakemake@log[['log']]
outplot1 <- snakemake@output[['plot1']]

theme <- snakemake@input[['theme']]
load(theme)
save.image(paste0(outplot1, '.Rdata'))

d <- read.table(ratios_f, sep="\t", header=T)
dnona <- d[!is.na(d$ratio),]
dnona$x <- ifelse(dnona$type =="PRIMARY", 'PRs', 'LMs')
dnona$x <- factor(dnona$x, levels=c('PRs', 'LMs'))


y_breaks <- guess_ticks(dnona$ratio)

p <- ggplot(data=dnona, aes(x=x, y=ratio))+geom_violin()+geom_boxplot(width=0.1, outlier.size=0.5)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
unmute_theme+ylab('SBS8/SBS1')

ggsave(outplot1, plot=p, width=89, height=89, units="mm")

wi <- wilcox.test(dnona[dnona$x=="LMs", 'ratio'], dnona[dnona$x=="PRs", 'ratio'], side="greater")
sink(log_f, append=TRUE)
print(wi)
print(wi$p.value)
sink()
