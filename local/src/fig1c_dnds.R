dnds_f  <- snakemake@input[['dnds']]
order_f  <- snakemake@input[['order']]
colors <- snakemake@input[['palette']]

#log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
library(ggpubr)
load(theme)

d <- read.table(dnds_f, sep="\t", quote="", header=FALSE, stringsAsFactors = TRUE)
colnames(d) <- c('model','estimate', 'upper', 'lower')
d$model <- gsub('_2nd', '', d$model, fixed=TRUE)


palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

#ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
#                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
#                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
#)
orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=TRUE)


d <- d[match(orderdf$model, d$model),]
if (!all(orderdf$model == d$model)) {
  stop('Issues in match between dnds data and model order in 1c')
}

d$model <- paste0(d$model, ifelse(!grepl('\\d$', d$model), '', ifelse(d$model=="CRC0282", 'PR', 'LM')))
names(pal) <- paste0(names(pal), ifelse(!grepl('\\d$', names(pal)), '', ifelse(names(pal)=="CRC0282", 'PR', 'LM')))
d$order <- seq(1, nrow(d))

y_breaks <- guess_ticks(d$lower,nticks=7)
y_breaks<- round(y_breaks, digits = 2)
print(y_breaks)
pdf('fig_1c_dnds.pdf')
p <- ggplot(d, aes(x=order(order, model), y=estimate, color=model)) +  geom_point(stat="identity", size=1.5) +
  geom_hline(yintercept=1,linetype=2,size=0.2)+
  geom_errorbar(aes(ymin=lower, ymax=upper, x=order(order,model)), width=0.3, size=0.3, color='black')+ylab('dN/dS estimate')+xlab('PDTs')+
  scale_color_manual(values=pal)+
  unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))# + ylim(NA,max(y_breaks))
print(p)
graphics.off()       
ggsave(outplot, plot=p, width=89, height=89, units="mm")
save.image(paste0(outplot, '.Rdata'))
