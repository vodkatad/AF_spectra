MR_f  <- snakemake@input[['mr']]
tree_f  <- snakemake@input[['tree']]
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
tree <- read.table(tree_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(tree) <- paste0(colnames(tree), "_tree")

m <- merge(mr, tree, by="row.names")
m$model <- m$Row.names

ci <- cor.test(m$mean_tree, m$mean_mr, method="spearman")
ci2 <- cor.test(m$mean_tree, m$mean_mr, method="pearson")


p <- ggplot(m, aes(x=mean_mr, y=mean_tree)) +  geom_point(aes(color=model), size=1) + geom_smooth(method='lm')+
  unmute_theme+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=pal)+xlab('MR')+ylab('SNV/indel trees dist')

ggsave(outplot, plot=p, width=89, height=56, units="mm")

sink(log_f)
print('n models')
print(nrow(m))
print('spearman:')
print(ci)
print('pearson:')
print(ci2)
sink()
save.image(paste0(outplot, '.Rdata'))

m <- m[m$model != "CRC0282PR",]

ci <- cor.test(m$mean_tree, m$mean_mr, method="spearman")
ci2 <- cor.test(m$mean_tree, m$mean_mr, method="pearson")
p <- ggplot(m, aes(x=mean_mr, y=mean_tree)) +  geom_point(aes(color=model), size=1) + geom_smooth(method='lm')+
  unmute_theme+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=pal)+xlab('MR')+ylab('SNV/indel trees dist')

sink(log_f, append=TRUE)
print('n models')
print(nrow(m))
print('spearman:')
print(ci)
print('pearson:')
print(ci2)
sink()

ggsave(paste0('noMSI_', outplot), plot=p, width=89, height=56, units="mm")
