library(ggplot2)
load('/scratch/trcanmed/AF_spectra/local/share/data/cn.Rdata')

nrow(our)==5*9 #???
our$sample <- rownames(our)
our <- our[,c(2,1)]
colnames(our) <- c('sample','dist_CNV') # subs all
#our$MR <- our$dist_CNV / 0.000000001
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
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

LEVEL <- 0.99
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'dist_CNV'], LEVEL) })
colnames(ic_clones) <- unique(our$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)
ggplot(pdata, aes(x=model, y=mean, color=model)) +  geom_point(position=position_dodge(1), stat="identity") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, size=1, position=position_dodge(1))+theme_bw()+ggtitle('CNV distance')+ylab('MEDICC2 events')

our$model_clone <- paste0(our$model, "_", our$clone)

#cbPalette <- c("#ff5733", "#ff7433","#ff8d33", 
#              #"#9d01fc","#a91efe", "#be52ff",
#              "#f607b9","#fb49ce","#f998e0",
#              "#155d00","#239203","#2fc603",
#              "#77a003","#95c805","#bcfc08",
#              "#0829fc","#4a62fb","#95a3fd"
#              )

#         282_01  282_05  282_07, 327_02   327_04  327_08   441_01 441_03 441_10  1078_2  1078_7 1078_9,     1307_02  1307_08 1307_08d 1307_09, 1502_3,1502_3A,1502_8,1502_9, 1502_9c, 1502_10, 1502_10b,1599LM1,1599lm3,1599pr1,1599pr10
colors = "#cc3300,#ff4000,#ff6633,#f607b9,#fb49ce,#f998e0,#9900ff,#ad33ff,#c266ff,#155d00,#239203,#2fc603,#77a003,#95c805,#95c805,#bcfc08,#0829fc,#0829fc,#4a62fb,#95a3fd,#95a3fd,#003399,#003399,#ff9900,#ffad33,#ffff00,#ffff66"


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

if (n == length(cbPalette)) {
  ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
    geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('CNV distance')+ylab('MEDICC2 events')+
    geom_point(data=our, aes(x=model, y=dist_CNV, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
    theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), legend.position="none", axis.title.y=element_text(size=15))+scale_color_manual(values=cbPalette)+scale_shape_manual(values=c(18,20))
} else {
  ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
    geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('CNV distance')+ylab('MEDICC2 events')+
    geom_point(data=our, aes(x=model, y=dist_CNV, color=model_clone, shape=time), stat="identity", size=3, position=position_dodge(0.2))+
    theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), legend.position="none", axis.title.y=element_text(size=15))+scale_shape_manual(values=c(18,20))
}
ggsave('cn_lengths.svg', height=5.25, width=5.25, units="in")

##
our2 <- read.table('/scratch/trcanmed/AF_spectra/dataset/vitro_gained_SNV', sep="\t", header=F, stringsAsFactors = F)
colnames(our2) <- c('sample','estimate')

our2$model <- sapply(our2$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})

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
ic_clones2 <- sapply(unique(our2$model), function(x) { confidence_interval(our2[our2$model==x,'estimate'], LEVEL) })
colnames(ic_clones2) <- unique(our2$model)
pdata2 <- as.data.frame(t(ic_clones2))
pdata2$model <- rownames(pdata2)


mm <- merge(pdata, pdata2, by="model")
m <- mm[,c('model','mean.x','mean.y')]
colors <- c("#cc3300,#f607b9,#9900ff,#155d00,#77a003,#0829fc,#ff9900,#ffff00")

name_x <- 'Distance_CN'
name_y <- 'Gained_SNV'
cbPalette <- unlist(strsplit(colors, ','))
colnames(m) <- c('model',name_x, name_y)
ci <- cor.test(m[, name_x], m[, name_y])

ggplot(m, aes_string(x=name_x, y=name_y)) +  geom_point(aes(color=model), size=3) + geom_smooth(method='lm')+
  theme_bw()+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=cbPalette)+theme(text = element_text(size = 15))
#ggplot(m, aes(x=MR_SNV, y=MR_indel)) +  geom_point(aes(color=model), size=3) + geom_smooth(method='lm', se=TRUE)+
#  theme_bw()+labs(caption=paste0('pearson=', ci$estimate, ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=cbPalette)

m <- m[m$model != 'CRC0282',]
cbPalette <- cbPalette[-1]
ci <- cor.test(m[, name_x], m[, name_y])
ggplot(m, aes_string(x=name_x, y=name_y)) +  geom_point(aes(color=model), size=3) + geom_smooth(method='lm')+
  theme_bw()+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=cbPalette)+theme(text = element_text(size = 15))

ggplot(m, aes_string(x=name_x, y=name_y)) +  geom_point(aes(color=model), size=3) + geom_smooth(method='lm')+
  theme_bw()+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=cbPalette)+theme(text = element_text(size = 15), legend.position="None")

ggsave('cn_cor.svg', height=5.25, width=5.25, units="in")
