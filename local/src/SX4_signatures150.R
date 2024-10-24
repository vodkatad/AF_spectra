library(ggplot2)
library(reshape)

outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
s30_f <- snakemake@input[['s30']]
s150_f <- snakemake@input[['s150']]

load(theme)

s150 <- read.table(s30_f, sep="\t", stringsAsFactors = F, header=T)
s30 <- read.table(s150_f, sep="\t", stringsAsFactors = F, header=T)

rownames(s150)[1:3] <- paste0(rownames(s150)[1:3], '-150x')
rownames(s30) <- paste0(rownames(s30), '-30x')
ss <- rbind(s150, s30)

# remove 282
ss <- ss[!grepl('282', rownames(ss)),]
colnames(ss) <- gsub('X', 'SBS', colnames(ss))
pd <- data.frame(SBS1=ss$SBS1, SBS8=ss$SBS8, SBS18=ss$SBS18)
pd$model <- substr(rownames(ss), 0, 7)
pd$class <- ifelse(grepl('bulk', rownames(ss)), 'Pre-existing', ifelse(grepl('30x', rownames(ss)), 'MA 30x', 'MA 150x'))

pdlong <- melt(data=pd, id=c("class", 'model'))
colnames(pdlong) <- c('class', 'model', 'Signature', 'exp')
pdlong$class <- factor(pdlong$class, c('Pre-existing', 'MA 30x', 'MA 150x'))

p <- ggplot(data=pdlong, aes(y=exp,fill=Signature, x=class))+
  geom_col(color='black')+
  facet_wrap(~model+Signature)+ylab('Relative contribution')+
  theme_bw() + theme(legend.position='bottom',
    text = element_text(size = textSize, family='sans'),
    axis.title = element_text(size = largerSize),
    axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
    axis.text.y = element_text(size = textSize, color="black"),
    plot.title = element_text(size = largerSize, hjust = 0.5),
    legend.title = element_text(size=largerSize, hjust = 0.5),
    legend.text = element_text(size=textSize),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(color = "black")
  )+
  xlab('')+scale_y_continuous(breaks=(seq(0, 0.6, by=0.15)),limits=c(-0.05, 0.6),expand = c(0, 0))

ggsave(outplot, plot=p, width=100*death_conversion_dpi96, height=100*death_conversion_dpi96, units="mm")