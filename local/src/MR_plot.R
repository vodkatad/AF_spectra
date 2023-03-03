#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
colors <- args[3]
ggtheme <- args[4]
load(ggtheme)

our <- read.table(infile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(our) <- c('sample','MR_edu')
our$MR <- our$MR_edu / 0.000000001
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((1+interval)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

LEVEL <- 0.99
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'MR'], LEVEL) })
colnames(ic_clones) <- unique(our$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)
ggplot(pdata, aes(x=model, y=mean, color=model)) +  geom_point(position=position_dodge(1), stat="identity") +
geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, size=1, position=position_dodge(1))+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')

our$model_clone <- paste0(our$model, "_", our$clone)

#cbPalette <- c("#ff5733", "#ff7433","#ff8d33", 
#              #"#9d01fc","#a91efe", "#be52ff",
#              "#f607b9","#fb49ce","#f998e0",
#              "#155d00","#239203","#2fc603",
#              "#77a003","#95c805","#bcfc08",
#              "#0829fc","#4a62fb","#95a3fd"
#              )
n <- length(levels(as.factor(our$model_clone)))
cbPalette <- unlist(strsplit(colors, ','))
#if (n == length(cbPalette)) {
# ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
# geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+
#   geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", shape=18, size=4, position=position_dodge(0.2))+
#   theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), legend.position="none", axis.title.y=element_text(size=15))+scale_color_manual(values=cbPalette)
# } else {
# ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
# geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+
#   geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", shape=18, size=4, position=position_dodge(0.2))+
#   theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), legend.position="none", axis.title.y=element_text(size=15))
# }

# ggsave(outfile)

# shape clones
our$time <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})


ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

if (n == length(cbPalette)) {
  ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
    geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
    geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
    ctheme+scale_color_manual(values=cbPalette)+scale_shape_manual(values=c(18,20))
} else {
  ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
    geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
    geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
    ctheme+scale_shape_manual(values=c(18,20))
}
ggsave(outfile)

save.image(paste0(outfile, '.Rdata'))

q('no')
### Reorderered svg for Andrea with theme
basename <- substr(outfile, 1, nchar(outfile)-4)
reordered <- paste0(basename, "_reordered.svg")
new_order <- c('CRC0282', 'CRC0327', 'CRC0441', 'CRC1502', 'CRC1599PR', 'CRC1599LM', 'CRC1078', 'CRC1307')
if (n == length(cbPalette)) {
  # We need to keep colors in track, we do this adding it to our before setting factors.
  ordered <- our$model_clone[order(our$model_clone)]
  palette_df <- data.frame(palette=cbPalette, model_clone=unique(ordered), stringsAsFactors=FALSE)# unique keeps the order, gasp
  #our <- merge(our, palette_df, by="model_clone")
  our$model <- factor(our$model, levels=new_order)
  pdata$model <- factor(pdata$model, levels=new_order)
  palette_values <- c(palette_df$palette)
  names(palette_values) <- palette_df$model_clone
  ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
    geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
    geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
    unmute_theme+ctheme+scale_color_manual(values=palette_values)+scale_shape_manual(values=c(18,20))
} else {
  ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
    geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
    geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", shape=18, size=4, position=position_dodge(0.2))+
    unmute_theme+ctheme+scale_shape_manual(values=c(18,20))
}
ggsave(reordered, height=5.25, width=5.25, units="in")


glibrary(ggsignif)


ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
  geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=cbPalette)+scale_shape_manual(values=c(18,20))+
  geom_signif(data=our, mapping=aes(x=model, y=MR), 
              comparisons = list(c("CRC0282", "CRC1307"),c('CRC1307','CRC1599PR')), map_signif_level=TRUE)

ggsave('MR_edu_SNV_ggsignif.svg', height=5.25, width=5.25, units="in")

ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
  geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=cbPalette)+scale_shape_manual(values=c(18,20))+
  geom_signif(data=our, mapping=aes(x=model, y=MR), 
              comparisons = list(c("CRC1599PR", "CRC1599LM"),c('CRC1078','CRC1307')), map_signif_level=TRUE)

ggsave('MR_edu_SNV_ggsignif2.svg', height=5.25, width=5.25, units="in")
