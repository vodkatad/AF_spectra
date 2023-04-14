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

tcga_p <- ggplot(muts_tcga, aes_string(x='status', y='burden', fill='status'))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(size=0.15)+
  scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+
  geom_signif(comparisons=list(c('MUT','WT')), size=0.2, textsize=1)+xlab('DNAH5')+ggtitle('TCGA')+
  unmute_theme

wes_p <- ggplot(muts_wes, aes_string(x='status', y='burden', fill='status'))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(size=0.15)+
  scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+
  geom_signif(comparisons=list(c('MUT','WT')), size=0.2, textsize=1)+xlab('DNAH5')+ggtitle('vWES')+
  unmute_theme

p <- ggarrange(tcga_p, wes_p,  common.legend = TRUE)


ggsave(outplot, plot=p, width=89, height=44, units="mm")
save.image(paste0(outplot, '.Rdata'))