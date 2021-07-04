datau <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1307_unpaired/platypus_nobin/all.MR_ov', sep="\t", header=T, stringsAsFactors = F)
datap <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1307_platypus_nobin/all.MR_ov', sep="\t", header=T, stringsAsFactors = F)
our <- datap[,c('end','MR_EDU')]

colnames(our) <- c('sample','MR_edu')
our$MR <- our$MR_edu / 0.000000001
our <- our[!grepl('-M', our$sample),]

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
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'MR'], LEVEL) })
colnames(ic_clones) <- unique(our$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)
pdatap <- pdata
ggplot(pdata, aes(x=model, y=mean, color=model)) +  geom_point(position=position_dodge(1), stat="identity") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, size=1, position=position_dodge(1))+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')

our$model_clone <- paste0(our$model, "_", our$clone)

ourp <- our

our <- datau[,c('end','MR_EDU')]

colnames(our) <- c('sample','MR_edu')
our <- our[!grepl('-M', our$sample),]

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
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
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
our$model <- paste0(our$model, '_unp')
pdata$model <- paste0(pdata$model, '_unp')
pdata <- rbind(pdatap, pdata)
our <- rbind(ourp, our)
our <- our[!grepl('-M', our$sample),]
ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+
  geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", size=2.5, position=position_dodge(0.2))+
  theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), axis.title.y=element_text(size=15))+scale_shape_manual(values=c(18,20))