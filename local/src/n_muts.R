library(ggplot2)

n_f  <- snakemake@input[['n']]
colors <- snakemake@input[['palette']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
outplotMSI <- snakemake@output[['plotMSI']]

#rule <- snakemake@rule

theme <- snakemake@input[['theme']]
load(theme)

save.image(paste0(outplot, '.Rdata'))

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model
names(pal) <- paste0(names(pal), ifelse(!grepl('\\d$', names(pal)), '', ifelse(names(pal)=="CRC0282", 'PR', 'LM')))

d <- read.table(n_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(d) <- c('#_SNVs', 'sample')

d$type <- ifelse(grepl('bulk', d$sample), 'parental', 'T0_clone')
d$model <- rep('', nrow(d))
d[d$type=="parental", 'model'] <- sapply(strsplit(d[d$type=="parental", 'sample'], "_"), '[[', 2)
d[d$type!="parental", 'model'] <- sapply(strsplit(d[d$type!="parental", 'sample'], "-"), '[[', 1)
#d$model <- ifelse(d$type=="parental", sapply(strsplit(d$sample, "_"), '[[', 2), sapply(strsplit(d$sample, "-"), '[[', 1))

d$color <- d$type
d <- d[order(d$model, d$type),]
j <- 0
oldmodel <- ''
for (i in seq(1, nrow(d))) {
  if (d[i, 'type'] ==  "T0_clone") {
    if (oldmodel == "" || oldmodel == d[i, 'model']) {
      j <- j + 1
    } else {
      j <- 1
    }
    d[i, 'type'] <- paste0('TO_clone_', j)
    oldmodel <- d[i, 'model']
  }
}

d$model <-  paste0(d$model, ifelse(!grepl('\\d$', d$model), '', ifelse(d$model=="CRC0282", 'PR', 'LM')))

dmsi <- d[d$model=="CRC0282PR",]
d <- d[d$model != "CRC0282PR",]

y_breaks <- guess_ticks(d$`#_SNVs`, fixed_max=24000)

d$x <- paste0(d$model,"_", d$type)
p <-ggplot(data=d, aes(fill=model, y=`#_SNVs`, x=x))+geom_col(position="dodge")+
      unmute_theme+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
    scale_fill_manual(values=pal)+
    scale_x_discrete(labels=d$color)+xlab('')

ggsave(outplot, plot=p, width=89, height=56, units="mm")

dmsi$x <- paste0(dmsi$model,"_", dmsi$type)
y_breaks <- guess_ticks(dmsi$`#_SNVs`, fixed_max=125000)

p <- ggplot(data=dmsi, aes(fill=model, y=`#_SNVs`, x=x))+geom_col(position="dodge")+
  unmute_theme+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  scale_fill_manual(values=pal)+
  scale_x_discrete(labels=d$color)+xlab('')

ggsave(outplotMSI, plot=p, width=89, height=56, units="mm")
