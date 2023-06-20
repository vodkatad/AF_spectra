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
guess_ticks_underzero <- function(values, nticks=5, fixed_max=NULL) {
  vmax <- max(values)
  if (is.null(fixed_max)) { 
    round_max <- vmax
  } else {
    round_max <- fixed_max
  }
  v_min<-min(values)
  my_breaks <- seq(v_min, round_max, length.out=nticks)
  return(my_breaks)
}
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
m$mean_mr<-m$mean_mr*1000000000
ci <- cor.test(m$mean_cn, m$mean_mr, method="spearman")
y_breaks <- guess_ticks(m$mean_cn+6)
x_breaks<-guess_ticks(m$mean_mr)
pdf('fig_4c_cor.pdf')
p <- ggplot(m, aes(x=mean_mr, y=mean_cn)) +  geom_point(aes(color=model), size=1) + geom_smooth(method='lm')+
  unmute_theme+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=pal)+xlab('MR')+ylab('MEDICC2 Events dist')+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)),expand=c(0,0))
print(p)
graphics.off()
ggsave(outplot, plot=p, width=89, height=56, units="mm")

sink(log_f)
print('n models')
print(nrow(m))
print(ci)
sink()
save.image(paste0(outplot, '.Rdata'))

m <- m[m$model != "CRC0282PR",]
y_breaks <- guess_ticks(m$mean_cn+7)

ci <- cor.test(m$mean_cn, m$mean_mr, method="spearman")

x_breaks<-guess_ticks(m$mean_mr)
print(x_breaks)

p <- ggplot(m, aes(x=mean_mr, y=mean_cn)) +  geom_point(aes(color=model), size=1) + geom_smooth(method='lm')+
  unmute_theme+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=pal)+xlab('MR')+ylab('MEDICC2 Events dist')+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0.01,max(x_breaks)),expand=c(0,0))

ggsave(paste0('noMSI_', outplot), plot=p, width=89, height=56, units="mm")
