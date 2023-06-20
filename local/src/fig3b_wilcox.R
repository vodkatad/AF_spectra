tcga_rds  <- snakemake@input[['TCGA']]
wes_rds  <- snakemake@input[['WES']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))
library(ggplot2)
library(ggpubr)

load(theme)

muts_tcga <- readRDS(tcga_rds)
muts_wes <- readRDS(wes_rds)

y_breaks <- guess_ticks(c(muts_tcga$burden, muts_wes$burden))
print(y_breaks)
pdf('fig_3b_wilcox_mb.pdf')
tcga_p <- ggplot(muts_tcga, aes_string(x='status', y='burden', fill='status'))+
  geom_boxplot(outlier.shape = NA,color='black')+geom_jitter(size=0.15, color="black", alpha=0.8)+
  stat_boxplot(geom ='errorbar', width = 0.3) +
  scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+
  #geom_signif(comparisons=list(c('MUT','WT')), size=0.2, textsize=1)
  xlab('TCGA')+ylab('Mutational burden')+labs(fill = "DNAH5")+
  unmute_theme+theme(axis.ticks.x = element_blank(),axis.text.x=element_blank()) +
  scale_y_continuous(breaks=y_breaks, limits=c(0, max(y_breaks)), expand = c(0, 0))
  
wes_p <- ggplot(muts_wes, aes_string(x='status', y='burden', fill='status'))+
  geom_boxplot(outlier.shape = NA,, color="black")+geom_jitter(size=0.15, color="black", alpha=0.8)+
  stat_boxplot(geom ='errorbar', width = 0.3) +
  scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+
  #geom_signif(comparisons=list(c('MUT','WT')), size=0.2, textsize=1)
  xlab('PDX')+ylab('Mutational burden')+labs(fill = "DNAH5")+
  unmute_theme+theme(axis.ticks.x = element_blank(),axis.text.x=element_blank()) +
  scale_y_continuous(breaks=y_breaks, limits=c(0, max(y_breaks)), expand = c(0, 0))

p <- ggarrange( wes_p,tcga_p,  common.legend = TRUE, legend="right") #,axis.text.x = element_blank(), axis.ticks.x = element_blank()

print(p)
graphics.off()


ggsave(outplot, plot=p,, width=89, height=89, units="mm") #
save.image(paste0(outplot, '.Rdata'))