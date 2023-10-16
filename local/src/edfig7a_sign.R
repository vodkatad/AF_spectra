library(reshape)
library(ggplot2)

signus_f  <- snakemake@input[['signus']]
signclevers_f  <- snakemake@input[['signclevers']]

log_f <- snakemake@log[['log']]
outplot1 <- snakemake@output[['plot1']]
outplot2 <- snakemake@output[['plot2']]

theme <- snakemake@input[['theme']]
load(theme)
save.image(paste0(outplot1, '.Rdata'))

wanted <- c('1', '8', '18')

data <- read.table(signus_f, sep="\t", header=TRUE)
wdata <- data[,colnames(data) %in% paste0('X', wanted)]
#normalized <- wdata/rowSums(wdata)
wdata$class <- sapply(strsplit(rownames(wdata), "_"), '[[', 2)
wdata$id <- rownames(wdata)
mp <- melt(wdata)
colnames(mp) <- c('class', 'id', 'Signature', 'exp')
mp$Signature <- gsub('X', 'SBS', mp$Signature, fixed=TRUE)
mp$class <- ifelse(mp$class=="bulk", 'Pre-existing', ifelse(mp$class=="vitroMA", "MA in vitro", "MA in vivo"))
mp$class <- factor(mp$class, c('Pre-existing', 'MA in vitro', 'MA in vivo'))
mp$Signature <- factor(mp$Signature, c('SBS1', 'SBS8', 'SBS18'))

p <- ggplot(data=mp, aes(y=exp,fill=Signature, x=class))+
    geom_boxplot(outlier.shape = NA,color='black')+
    geom_jitter(size=0.15, color="black", alpha=0.8)+
    facet_wrap(~Signature)+ylab('Relative contribution')+
    xlab('')+unmute_theme+scale_y_continuous(breaks=(seq(0, 0.5, by=0.1)))+theme(legend.position="bottom")

death_conversion_dpi96 = 96/72
ggsave(outplot1, plot=p, width=78.5*death_conversion_dpi96, height=71.9*death_conversion_dpi96, units="mm")

data <- read.table(signclevers_f, sep="\t", header=TRUE)
wdata <- data[,colnames(data) %in% paste0('X', wanted)]
#normalized <- wdata/rowSums(wdata)
wdata$class <- sapply(strsplit(rownames(wdata), "_"), '[[', 2)
wdata$id <- rownames(wdata)
mp <- melt(wdata)
colnames(mp) <- c('class', 'id', 'Signature', 'exp')
mp$Signature <- gsub('X', 'SBS', mp$Signature, fixed=TRUE)
mp$class <- ifelse(mp$class=="leaves", 'Leaves', 'Truncal')
mp$class <- factor(mp$class, c('Truncal', 'Leaves'))
mp$Signature <- factor(mp$Signature, c('SBS1', 'SBS8', 'SBS18'))

p <- ggplot(data=mp, aes(y=exp,fill=Signature, x=class))+
    geom_boxplot(outlier.shape = NA,color='black')+
    geom_jitter(size=0.15, color="black", alpha=0.8)+
    facet_wrap(~Signature)+ylab('Relative contribution')+
    xlab('')+unmute_theme+scale_y_continuous(breaks=(seq(0, 0.5, by=0.1)))+theme(legend.position="bottom")


ggsave(outplot2, plot=p, width=78.5*death_conversion_dpi96, height=71.9*death_conversion_dpi96, units="mm")

# Further work
# remove ticks
# remove upper grey boxes?