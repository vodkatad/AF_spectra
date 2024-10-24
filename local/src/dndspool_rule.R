library(ggplot2)

order_f  <- snakemake@input[['order']]
colors <- snakemake@input[['palette']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
dnds_subcl_f <- snakemake@input[['dnds_subcl']]
dnds_pool_f <- snakemake@input[['dnds_pool']]

#order_f  <- '/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/sorted_MR_avg.tsv'
#colors <- '/scratch/trcanmed/AF_spectra/local/share/data/model_palette.rds'
#outplot <- snakemake@output[['plot']]
#theme <- '/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_5.Rdata'
load(theme)

palette_df <- readRDS(colors)
pal <- c(palette_df$palette, 'darkgrey')
names(pal) <- c(palette_df$model, 'pooled_MSS')

dndsall <- read.table(dnds_subcl_f, sep="\t", header=FALSE)
colnames(dndsall) <- c('model','estimate', 'lower', 'upper')

dndspooled <- read.table(dnds_pool_f, sep="\t", header=TRUE)
colnames(dndspooled) <- c('model','estimate', 'lower', 'upper')

pooldf <- cbind('pooled_MSS', dndspooled['wall', c('estimate', 'lower', 'upper')])
colnames(pooldf)[1] <- 'model'
              

d <- rbind(dndsall, pooldf)


orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=TRUE)
orderdf <- rbind(orderdf, data.frame(lower=1, upper=1, mean=1, model='pooled_MSS'))


d <- d[match(orderdf$model, d$model),]
if (!all(orderdf$model == d$model)) {
  stop('Issues in match between dnds data and model order in 1c')
}

d$model <- paste0(d$model, ifelse(!grepl('\\d$', d$model), '', ifelse(d$model=="CRC0282", 'PR', 'LM')))
names(pal) <- paste0(names(pal), ifelse(!grepl('\\d$', names(pal)), '', ifelse(names(pal)=="CRC0282", 'PR', 'LM')))
d$order <- seq(1, nrow(d))

y_breaks <- guess_ticks(c(d$lower, d$upper, d$estimate),nticks=7)
y_breaks<- round(y_breaks, digits = 2)

# p <- ggplot(d, aes(x=order(order, model), y=estimate, color=model)) +  geom_point(stat="identity", size=1.5) +
#   geom_hline(yintercept=1,linetype=2,size=0.2)+
#   geom_errorbar(aes(ymin=lower, ymax=upper, x=order(order,model)), width=0.3, size=0.3, color='black')+ylab('dN/dS estimate')+xlab('PDTs')+
#   scale_color_manual(values=pal)+
#   unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
#                      axis.ticks.x = element_blank(),
#                      legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+
#   scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))# + ylim(NA,max(y_breaks))
# 
# print(p)

p <- ggplot(d, aes(x=order(order, model), y=estimate, color=model)) +  geom_point(stat="identity") +
  geom_hline(yintercept=1,linetype=2,size=0.2)+
  geom_errorbar(aes(ymin=lower, ymax=upper, x=order(order,model)), width=0.3, size=0.3, color='black')+ylab('dN/dS estimate')+xlab('PDTs')+
  scale_color_manual(values=pal)+
  unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))# + ylim(NA,max(y_breaks))


ggsave(outplot, plot=p, width=89, height=89, units="mm")