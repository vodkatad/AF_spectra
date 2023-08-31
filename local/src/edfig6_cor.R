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
m$mean_tree<- m$mean_tree
m$mean_mr<-m$mean_mr*1000000000
y_breaks <- guess_ticks(m$mean_tree,fixed_max=120000)
x_breaks<- guess_ticks(m$mean_mr)

p <- ggplot(m, aes(x=mean_mr, y=mean_tree)) + geom_smooth(method='lm') +geom_point(aes(color=model), size=1) +
  scale_y_continuous(breaks=y_breaks, limits=c(-15000,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)),expand=c(0,0))+
  unmute_theme + scale_color_manual(values=pal)+xlab('MR [SNV/(Gbp*division)]')+ylab('SNV')+theme(legend.position="none",legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))#+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4)))

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
  unmute_theme + scale_color_manual(values=pal)+xlab('MR [SNVs/(Gbp*division)]')+ylab('MEDICC2 events')#caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))

sink(log_f, append=TRUE)
print('n models')
print(nrow(m))
print('spearman:')
print(ci)
print('pearson:')
print(ci2)
sink()

ggsave(paste0('noMSI_', outplot), plot=p, width=89, height=56, units="mm")
