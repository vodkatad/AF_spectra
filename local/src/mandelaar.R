library(ggplot2)

sourcedata_f  <- snakemake@input[['sourcedata']]
log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]

theme <- snakemake@input[['theme']]
load(theme)
save.image(paste0(outplot, '.Rdata'))


d <- read.table(sourcedata_f, sep="\t", header=T)
d$type <- ifelse(grepl('1$', d$ori), 'Mets', 'PRs')
pd <- d[grepl('SBS1', d$ori),]
ggplot(data=pd, aes(x=totals, color=type))+geom_density()
pd8 <- d[grepl('SBS8', d$ori),]
ggplot(data=pd8, aes(x=totals, color=type))+geom_density()

# Here we assume that patients are listed in the same order for both signatures in the source data
# I found no way to prove this but it makes sense...
pd8$ratio <- pd8$totals / pd$totals
pd8$type <- factor(pd8$type, levels=c('PRs', 'Mets'))
pd8orig <- pd8
pd8 <- pd8[!is.na(pd8$ratio),]
y_breaks <- guess_ticks(pd8$ratio)


p <- ggplot(data=pd8, aes(y=ratio, x=type))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0, size=0.15)+
  stat_boxplot(geom ='errorbar', width = 0.3) +
  xlab('')+ylab('ratio')+
  unmute_theme+theme(axis.ticks.x = element_blank()) +#,axis.text.x=element_blank()
  scale_y_continuous(breaks=y_breaks, limits=c(0, max(y_breaks)), expand = c(0, 0))

wilc <- wilcox.test(pd8[pd8$type=='Mets', 'ratio'], pd8[pd8$type=='PRs', 'ratio'], alt="greater")

sink(log_f)
print('Total=')
print(table(pd8$type))
print('Removed SBS1==0, two times')
print(sum(pd$totals == 0))
print(sum(is.na(pd8orig$ratio)))
print(wilc)
print(wilc$p.value)
sink()

ggsave(outplot, plot=p, width=89, height=89, units="mm")


