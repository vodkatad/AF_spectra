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

tcga_p <- ggplot(muts_tcga, aes_string(x='status', y='burden', fill='status'))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(size=0.15, color="gray34", alpha=0.8)+
  scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+
  geom_signif(comparisons=list(c('MUT','WT')), size=0.2, textsize=1)+xlab('DNAH5')+ggtitle('TCGA')+
  unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(0, max(y_breaks)), expand = c(0, 0))
  
wes_p <- ggplot(muts_wes, aes_string(x='status', y='burden', fill='status'))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(size=0.15, color="gray34", alpha=0.8)+
  scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+
  geom_signif(comparisons=list(c('MUT','WT')), size=0.2, textsize=1)+xlab('DNAH5')+ggtitle('PDX WES')+
  unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(0, max(y_breaks)), expand = c(0, 0))

p <- ggarrange(tcga_p, wes_p,  common.legend = TRUE, legend="right")


ggsave(outplot, plot=p, width=89, height=89, units="mm")
save.image(paste0(outplot, '.Rdata'))