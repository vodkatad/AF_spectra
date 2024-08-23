#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
mutfile <- args[1]
theme <- args[2]
maxcn <- as.numeric(args[3])
sample_name<-strsplit(mutfile,split='.',fixed=TRUE)[[1]][1]

outputplot <- args[4]
outputobj <- args[5]

load(theme)

muts <- read.table(gzfile(mutfile), header=FALSE, stringsAsFactors=FALSE)
colnames(muts) <- c("chr","b","e","id")
muts$af <- as.numeric(sapply(muts$id, function(x) { strsplit(x, ":")[[1]][7]}))
muts$cn <- sapply(muts$id, function(x) { strsplit(x, ":")[[1]][8]})
muts$cn <- as.numeric(muts$cn)
muts <- muts[muts$cn <= maxcn,]#
#due<-c(1/2)
#tre<-c(1/3,2/3)
#quattro<-c(1/4,1/2,3/4)
#cinque<-c(1/5,2/5,3/5,4/5)
#muts$hline<-
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
  geom_histogram(aes(fill=cn), alpha=0.3, binwidth=0.05, position = 'identity')+
  
  guides(fill =guide_legend(title='CN'),color = guide_legend(title = "CN")) 
  #geom_histogram(fill='white', alpha=0.3, binwidth=0.05, position = 'identity')#+
  #geom_density(position="identity", aes(y=after_stat(count*0.05)), bw = 0.05)
  #scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  #scale_x_continuous(breaks=x_breaks,limits=c(0, max(x_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  #unmute_theme#+theme(legend.position="none", axis.text.x = element_blank(), 
        	 #      axis.ticks.x = element_blank(),
                 #     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))      

ggp <- ggplot_build(p)
#maxy <- max(ggp$data[[2]]$count * 0.05) # [[2]] if histogram is removed
maxy <- max(ggp$data[[1]]$y+5) # [[1]] is histogram is removed
y_breaks <- guess_ticks(maxy)
x_breaks <- guess_ticks(muts$af, fixed_max=1)
y_labels<-c('0',as.character(round(y_breaks[2])),as.character(round(y_breaks[3])),as.character(round(y_breaks[4])),as.character(round(y_breaks[5])))
vline_list <- list(
  '2' = c(1/2),
  '3' = c(1/3,2/3),
  '4' = c(1/4,1/2,3/4),
  '5'=c(1/5,2/5,3/5,4/5)
)
vline_data <- do.call(rbind, lapply(names(vline_list), function(cn) {
  data.frame(cn = cn, vline = vline_list[[cn]])
}))

p <- p + scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0),labels=y_labels)+# + ylim(min(y_breaks),max(y_breaks))+
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),limits=c(0, 1),expand = c(0, 0),labels=c('0','0.25','0.5','0.75','1'))+
  geom_vline(data = vline_data, aes(xintercept = vline), color = "red", linetype = "dashed")+
  facet_grid(. ~ cn)+# + ylim(min(y_breaks),max(y_breaks))+ #
  ylab(sample_name)+xlab('VAF')+
  
  theme(strip.background = element_blank(),strip.text.x = element_blank())+  theme(legend.key.height= unit(0.2, 'cm'),
        legend.key.width= unit(0.2, 'cm'))+theme(legend.title=element_blank())+

  unmute_theme#+theme(legend.position="none", axis.text.x = element_blank(), 
       	 #      axis.ticks.x = element_blank(),
                 #     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))      

save.image('p.Rdata')
ggsave(file=outputplot, plot=p,width = 180,height=36,units = "mm",dpi = 300)

saveRDS(p, file=outputobj)