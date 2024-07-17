#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
mutfile <- args[1]
theme <- args[2]
maxcn <- as.numeric(args[3])
outputplot <- args[4]

load(theme)

muts <- read.table(gzfile(mutfile), header=FALSE, stringsAsFactors=FALSE)
colnames(muts) <- c("chr","b","e","id")
muts$af <- as.numeric(sapply(muts$id, function(x) { strsplit(x, ":")[[1]][7]}))
muts$cn <- sapply(muts$id, function(x) { strsplit(x, ":")[[1]][8]})
muts$cn <- as.numeric(muts$cn)
muts <- muts[muts$cn <= maxcn,]
muts$cn <- factor(muts$cn, levels=seq(1, maxcn))

max_y <- 0
all_cn <- levels(muts$cn)
find_max <- function(x, data) {
	myd <- data[data$cn==x,]
	return(nrow(myd))
}

#maxnvar <- sapply(all_cn, find_max, muts)

## TODO: get reasonable max studying after_stat output for density
# remove histogram
# put together with pdftk:             pdftk {params.plots}/*pdf cat output {output}
# https://stackoverflow.com/questions/68088769/overlay-kde-and-filled-histogram-with-ggplot2-r
# kde='true' seaborn plothisto
# https://stackoverflow.com/questions/11404531/r-ggplot2-adding-count-labels-to-histogram-with-density-overlay

p <- ggplot(data=muts, aes(x=af, color=cn)) +
  #geom_histogram(aes(y=after_stat(count), fill=cn), alpha=0.3, binwidth=0.05, position = 'identity')+
  geom_density(position="identity", aes(y=after_stat(count*0.05)), bw = 0.05)
  #scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  #scale_x_continuous(breaks=x_breaks,limits=c(0, max(x_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  #unmute_theme#+theme(legend.position="none", axis.text.x = element_blank(), 
        	 #      axis.ticks.x = element_blank(),
                 #     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))      

ggp <- ggplot_build(p)
maxy <- max(ggp$data[[1]]$count * 0.05) # [[2]] if histogram is removed
y_breaks <- guess_ticks(maxy)
x_breaks <- guess_ticks(muts$af, fixed_max=1)
#p <- p + scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
#  scale_x_continuous(breaks=x_breaks,limits=c(0, max(x_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
#  unmute_theme#+theme(legend.position="none", axis.text.x = element_blank(), 
        	 #      axis.ticks.x = element_blank(),
                 #     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))      

save.image('p.Rdata')
ggsave(file=outputplot, plot=p)

