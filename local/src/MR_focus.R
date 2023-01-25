
load('/scratch/trcanmed/AF_spectra/dataset/MR_edu_SNV.png.Rdata')

w <- c('CRC1078','CRC1307','CRC1599LM','CRC1599PR')
pdata <- pdata[pdata$model %in% w,]
grepps <- unlist(sapply(w, function(x) { our$model[grepl(x, our$model,)] }))
our <- our[our$model %in% grepps,]

cbPalette <- c('#155d00','#239203','#2fc603','#77a003','#95c805','#95c805','#bcfc08','#bcfc08','#ffff00','#ffff66','#ff9900','#ffad33')

pdata$model <- factor(pdata$model, levels=c('CRC1078','CRC1307','CRC1599PR','CRC1599LM'))

ggsave('MR_focus_time_size.png', width=7.574, height=7.574, units='in')

library(ggsignif)

ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15), axis.title.x = element_text(size=20),
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
  geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=cbPalette)+scale_shape_manual(values=c(18,20))+geom_signif(data=our, aes(x=model, y=MR),
    comparisons = list(c("CRC1078", "CRC1307"), c('CRC1599PR','CRC1599LM')),
    map_signif_level = TRUE, 
    #test="t.test"
  )

ggsave('MR_focus_time_size.png', width=7.574, height=7.574, units='in')


wilcox.test(our[grepl('CRC1078',our$sample),'MR_edu'], our[grepl('CRC1307',our$sample),'MR_edu'])

wilcox.test(our[grepl('CRC1599LM',our$sample),'MR_edu'], our[grepl('CRC1599PR',our$sample),'MR_edu'])


load('/scratch/trcanmed/AF_spectra/dataset/MR_edu_SNV.png.Rdata')

ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15), axis.title.x = element_text(size=20),
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

w <- c('CRC1502')
pdata <- pdata[pdata$model %in% w,]
grepps <- unlist(sapply(w, function(x) { our$model[grepl(x, our$model,)] }))
our <- our[our$model %in% grepps,]
cbPalette <- c('#0829fc','#0829fc','#4a62fb','#4a62fb','#95a3fd','#95a3fd','#003399','#003399')

ggplot(our, aes(x=time, y=MR, shape=time, color=model_clone)) +  geom_point(stat="identity", size=3, position=position_dodge(0.2)) +
  ctheme+scale_shape_manual(values=c(18,20))+scale_color_manual(values=cbPalette)+geom_signif(color='black',comparisons = list(c("1", "2")),
                                                                                              map_signif_level = TRUE, 
                                                                                              #test="t.test"
  )

ggsave('evil_clones_size.png', width=7.574, height=7.574, units='in')


wilcox.test(our[grepl('1',our$time),'MR_edu'], our[grepl('2',our$time),'MR_edu'])

#### mut subcl burden

our <- read.table('/scratch/trcanmed/AF_spectra/dataset/subclonal_SNV', sep="\t", header=F, stringsAsFactors = F)
colnames(our) <- c('sample','subclonal')
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)

our$subclonal <- our$subclonal * 1000000

our <- our[our$model !="CRC1307LMO", ]#should not be here!
our <- our[!grepl('-M', our$sample, fixed=T),]

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

LEVEL <- 0.99
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'subclonal'], LEVEL) })
colnames(ic_clones) <- unique(our$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

w <- c('CRC1078','CRC1307','CRC1599LM','CRC1599PR')
pdata <- pdata[pdata$model %in% w,]
grepps <- unlist(sapply(w, function(x) { our$model[grepl(x, our$model,)] }))
our <- our[our$model %in% grepps,]
our$time <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

                #1078_2  1078_7      1078_9   1307_02   1307_08   1307_09   1599LM1,1599lm3,  1599lm7 1599pr1,1599pr10, 
cbPalette <- c('#155d00','#239203','#2fc603','#77a003','#95c805','#bcfc08','#ff9900','#ffad33','#ff9900','#ffff00','#ffff66')
#cbPalette <- c('#ff9900','#ffad33','#ffff99', '#ffff00','#ffff66','#155d00','#239203','#2fc603','#77a003','#95c805','#bcfc08')



new_order <- c('CRC1599PR','CRC1599LM','CRC1078','CRC1307')
#our$model <- factor(our$model, levels=c('CRC1599PR','CRC1599LM','CRC1078','CRC1307'))
#our <- our[order(our$model),]

ordered <- our$model_clone[order(our$model_clone)]
palette_df <- data.frame(palette=cbPalette, model_clone=unique(ordered), stringsAsFactors=FALSE)# unique keeps the order, gasp
#our <- merge(our, palette_df, by="model_clone")
our$model <- factor(our$model, levels=new_order)
pdata$model <- factor(pdata$model, levels=new_order)
palette_values <- c(palette_df$palette)
names(palette_values) <- palette_df$model_clone


ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15), axis.title.x = element_text(size=20),
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

# att time, right now we have 0-1 instead of 0-1-2 TODO FIXME
ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('Subclonal muts')+ylab('mut/Mbp')+xlab('')+
  geom_point(data=our, aes(x=model, y=subclonal, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=palette_values)+scale_shape_manual(values=c(18,20))+geom_signif(data=our, aes(x=model, y=subclonal),
                                                                                              comparisons = list(c("CRC1078", "CRC1307"), c('CRC1599PR','CRC1599LM')),
                                                                                              map_signif_level = TRUE
                                                                                              #test="t.test"
  )

ggsave('subclonal_manual_size_correctcolors.pdf', width=7.574, height=7.574, units='in')
