dnds_f  <- snakemake@input[['dnds']]
colors <- snakemake@input[['palette']]

#log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
library(ggpubr)
load(theme)

d <- read.table(dnds_f, sep="\t", quote="", header=FALSE, stringsAsFactors = FALSE)
colnames(d) <- c('tmodel','estimate', 'lower', 'upper') 

d$model <- sapply(d$tmodel, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})
d$treat <- sapply(d$tmodel, function(x) {y<-strsplit(x, '_')[[1]][1]; return(y[1])})
d <- d[order(d$model, d$treat),]
d$order <- seq(1, nrow(d))

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

y_breaks <- guess_ticks(c(d$estimate, d$upper, d$lower), nticks=7)
y_breaks<- round(y_breaks, digits = 2)

p <- ggplot(d, aes(x=factor(order), y=estimate, color=model)) +  geom_point(stat="identity", size=1.5) +
  geom_hline(yintercept=1,linetype=2,size=0.2)+
  geom_errorbar(aes(ymin=lower, ymax=upper, x=order(order,model)), width=0.3, size=0.3, color='black')+ylab('dN/dS estimate')+xlab('')+
  scale_color_manual(values=pal)+
  unmute_theme+theme(legend.position="none",
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand = c(0, 0))+
  scale_x_discrete(breaks=c(1,2,3,4), labels=c('1'="NT", '2'="Afatinib", '3'="NT", '4'="Afatinib"))

ggsave(outplot, plot=p, width=89, height=89, units="mm")
save.image(paste0(outplot, '.Rdata'))
