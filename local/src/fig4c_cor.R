MR_f  <- snakemake@input[['mr']]
CN_f  <- snakemake@input[['cn']]
colors <- snakemake@input[['palette']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
data_f <- snakemake@output[['avgdata']]

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model
names(pal) <- paste0(names(pal), ifelse(!grepl('\\d$', names(pal)), '', ifelse(names(pal)=="CRC0282", 'PR', 'LM')))

mr <- read.table(MR_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mr) <- paste0(colnames(mr), "_mr")
cn <- read.table(CN_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(cn) <- paste0(colnames(cn), "_cn")

m <- merge(mr, cn, by="row.names")
m$model <- m$Row.names

ci <- cor.test(m$mean_cn, m$mean_mr, method="spearman")


p <- ggplot(m, aes(x=mean_mr, y=mean_cn)) +  geom_point(aes(color=model), size=1) + geom_smooth(method='lm')+
  unmute_theme+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=pal)+xlab('MR')+ylab('MEDICC2 Events dist')

ggsave(outplot, plot=p, width=89, height=56, units="mm")

sink(log_f)
print('n models')
print(nrow(m))
print(ci)
sink()
save.image(paste0(outplot, '.Rdata'))

m <- m[m$model != "CRC0282PR",]

ci <- cor.test(m$mean_cn, m$mean_mr, method="spearman")
p <- ggplot(m, aes(x=mean_mr, y=mean_cn)) +  geom_point(aes(color=model), size=1) + geom_smooth(method='lm')+
  unmute_theme+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=pal)+xlab('MR')+ylab('MEDICC2 Events dist')

ggsave(paste0('noMSI_', outplot), plot=p, width=89, height=56, units="mm")
