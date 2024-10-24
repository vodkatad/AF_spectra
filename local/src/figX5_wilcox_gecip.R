mb_crc_f  <- snakemake@input[['mb_crc']]
mb_panc_f  <- snakemake@input[['mb_panc']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
log_f <- snakemake@log[['log']]

save.image(paste0(outplot, '.Rdata'))
library(ggplot2)
library(ggpubr)

load(theme)

mb_crc <- read.table(mb_crc_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
mb_panc <- read.table(mb_panc_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mb_crc) <- c('status', 'burden')
colnames(mb_panc) <- c('status', 'burden')
mb_crc$status <- toupper(mb_crc$status)
mb_panc$status <- toupper(mb_panc$status)

o_wicrc <- wilcox.test(mb_crc[mb_crc$status=='MUT', 'burden'], mb_crc[mb_crc$status=='WT', 'burden'], alternative = 'greater')
o_wipanc <- wilcox.test(mb_panc[mb_panc$status=='MUT', 'burden'], mb_panc[mb_panc$status=='WT', 'burden'], alternative = 'greater')

mb_crc$burden <- log(mb_crc$burden + 0.01) # manually checked 0.03 was the min != 0
mb_panc$burden <- log(mb_panc$burden + 0.01) # manually checked 0.03 was the min != 0
y_breaks <- guess_ticks(c(mb_crc$burden, mb_panc$burden), fixed_min=-5) # again manual check
crc_p <- ggplot(mb_crc, aes(x=factor(status,level=c('WT','MUT')), y=burden))+#, fill='status'
  geom_boxplot(outlier.shape = NA,color='black')+#geom_jitter(size=0.15, color="black", alpha=0.8, height=0)+
  stat_boxplot(geom ='errorbar', width = 0.3) +
  xlab('CRC')+ylab('Mutational burden')+
  unmute_theme+theme(axis.ticks.x = element_blank()) +#,axis.text.x=element_blank()
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks), max(y_breaks)), expand = c(0, 0))
  
panc_p <- ggplot(mb_panc,aes(x=factor(status,level=c('WT','MUT')), y=burden))+#, fill='status'
  geom_boxplot(outlier.shape = NA, color="black")+#geom_jitter(size=0.15, color="black", alpha=0.8, height=0)+
  stat_boxplot(geom ='errorbar', width = 0.3) +
  xlab('Pan-cancer')+ylab('Mutational burden')+labs(fill = "DNAH5")+
  unmute_theme+theme(axis.ticks.x = element_blank())+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks), max(y_breaks)), expand = c(0, 0))

p <- ggarrange( crc_p,panc_p,  common.legend = TRUE, legend="none") #,axis.text.x = element_blank(), axis.ticks.x = element_blank()


wicrc <- wilcox.test(mb_crc[mb_crc$status=='MUT', 'burden'], mb_crc[mb_crc$status=='WT', 'burden'], alternative = 'greater')
wipanc <- wilcox.test(mb_panc[mb_panc$status=='MUT', 'burden'], mb_panc[mb_panc$status=='WT', 'burden'], alternative = 'greater')

sink(log_f)
print('CRC')
table(mb_crc$status)
print('panc')
table(mb_panc$status)
print('Wilcox on log')
print(wicrc)
print(wicrc$p.value)
print(wipanc)
print(wipanc$p.value)
print('wilcox on orig') # ovviamente uguali
print(o_wicrc)
print(o_wicrc$p.value)
print(o_wipanc)
print(o_wipanc$p.value)
sink()

ggsave(outplot, plot=p, width=178, height=89, units="mm") #
save.image(paste0(outplot, '.Rdata'))