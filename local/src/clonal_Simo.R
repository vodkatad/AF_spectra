setwd('/scratch/trcanmed/AF_spectra/datasetV2')

palette_df <- readRDS("/scratch/trcanmed/AF_spectra/local/share/data/palette.rds")
pal <- palette_df$palette
names(pal) <- palette_df$model_clone

palette_df <- readRDS("/scratch/trcanmed/AF_spectra/local/share/data/model_palette.rds")
pal2 <- palette_df$palette
names(pal2) <- palette_df$model_clone

subclonal <- read.table('vitro_gained_binSNV', header=FALSE, sep="\t", stringsAsFactors=FALSE)
#subclonal <- read.table('MR_edu_binSNV', header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(subclonal) <- c('sample', 'MR')
subclonal$MR_edu_cl <- subclonal$MR 
#subclonal$MR_edu_cl <- subclonal$MR / 0.000000001
MR <- read.table('vitro_gained_SNV', header=FALSE, sep="\t", stringsAsFactors=FALSE)
#MR <- read.table('MR_edu_SNV', header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(MR) <- c('sample', 'MR')
MR$MR_edu <- MR$MR 
#MR$MR_edu <- MR$MR / 0.000000001

merged <- merge(subclonal, MR, by='sample')
merged$smodel <- substr(merged$sample, 0, 7)
wanted <- c('CRC1599', 'CRC1307', 'CRC1078')
merged <- merged[merged$smodel != 'CRC0282',]

merged$model <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
merged$PDT <- paste0(merged$model, ifelse(grepl('\\d$', merged$model), ifelse(grepl('282', merged$model), 'PR', 'LM'), ''))
#merged$PDT <- factor(merged$PDT, levels=c('CRC1078LM','CRC1307LM','CRC1599PR','CRC1599LM'))

merged$clone <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
merged$clone2 <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
merged$model_clone <- paste0(merged$model, "_", merged$clone)
merged$time <- sapply(merged$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

# p <- ggplot(data=merged, aes(x=MR_edu, y=n)) +
#   geom_point(aes(color=model_clone, fill=model_clone), stat="identity", size=1, position=position_dodge(0.2))+
#   unmute_theme+scale_color_manual(values=pal, guide="none")+
#   scale_fill_manual(values=pal, guide="none")+
#   xlab('MR')+ylab('# subclonal SNVs')+theme(legend.position="right")

# average errors ############
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- sd(vector)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

LEVEL <- 0.99
ic_clones <- sapply(unique(merged$model), function(x) { confidence_interval(merged[merged$model==x,'MR_edu'], LEVEL) })
colnames(ic_clones) <- unique(merged$model)
pdataMR <- as.data.frame(t(ic_clones))
pdataMR$model <- rownames(pdataMR)


ic_clones <- sapply(unique(merged$model), function(x) { confidence_interval(merged[merged$model==x,'MR_edu_cl'], LEVEL) })
colnames(ic_clones) <- unique(merged$model)
pdatasub <- as.data.frame(t(ic_clones))
pdatasub$model <- rownames(pdatasub)

pdata <- merge(pdataMR, pdatasub, by='model')
# 
ggplot() +
  geom_point(data=merged, aes(x=MR_edu, y=MR_edu_cl, color=model_clone), stat="identity", size=2)+
  theme_bw()+
  geom_point(data=pdata, aes(x=mean.x, y=mean.y, shape=model, fill=model), size=2)+
  geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, ymin=lower.y, ymax=upper.y), width=.05, size=.2)+
  geom_errorbar(data=pdata, aes(x=mean.x, y=mean.y, xmin=lower.x, xmax=upper.x), width=.05, size=.2)+
  scale_fill_manual(values=pal2)+scale_color_manual(values=pal, guide='none')+
  #xlab('# SNVs gained')+ylab('# clonal SNVs gained')+theme(legend.position="right")+geom_abline(slope=1, intercept=0)
  xlab('MR')+ylab('clonal MR')+theme(legend.position="right")+geom_abline(slope=1, intercept=0)


### n bulk muts on CN1-2-3
d <- read.table(gzfile('/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/tree/bulk.var_cnv.tsv.gz'), sep="\t", header=F, stringsAsFactors = F)
colnames(d) <- c('chr', 'b', 'e', 'info')

infos <- strsplit(d$info, split=':', fixed=T)
d$af <- as.numeric(sapply(infos, '[', 7))
d$cn <- as.numeric(sapply(infos, '[', 8))
d$ref <- as.numeric(sapply(infos, '[', 5))
d$alt <- as.numeric(sapply(infos, '[', 6))
d$tot <- d$ref+d$alt
summary(d$tot)
d$vaf <- d$alt / d$tot
d <- d[d$cn %in% c(1,2,3),]
ggplot(data=d, aes(x=vaf, y=..count.., color=as.factor(cn)))+geom_density(position="identity")+theme_bw()
d$caf <- d$vaf * d$cn
ggplot(data=d, aes(x=caf, y=..count..))+geom_density(position="identity")+theme_bw()

                                                           

rbinom <- function(mut) {
  refreads <- as.integer(mut[7])
  mutreads <- as.integer(mut[8])
  cn <- as.integer(mut[6])
  tot <- refreads+mutreads
  s <- 1/cn
  sqrt_1_s <- (1-s)**(1/2)
  sqrt_s <- (s)**(1/2)
  sqrt_n <- (tot)**(1/2)
  if (cn==0) {
    return(0)
  }
  #print(paste0('Cn= ', cn))
  #print(paste0('refreads= ', refreads))
  #print(paste0('mutreads= ', mutreads))
  #s - (sqrt(1-s)*sqrt(s))/(sqrt(n)  < x < s + (sqrt(1-s)*sqrt(s))/(sqrt(n)
  if (mutreads/tot >  s - ((sqrt_1_s*sqrt_s)/sqrt_n) && mutreads/tot <  s + ((sqrt_1_s*sqrt_s)/sqrt_n)) {
    return(1)
  } else {
    return(0)
  }
}
new_binom <- apply(d, 1, rbinom)
d$bino <- new_binom

db <- d[d$bino ==1,]
dnb <- d[d$bino ==0,]
ggplot(data=db, aes(x=caf, y=..count..))+geom_density(position="identity")+theme_bw()

ggplot(data=d, aes(x=caf))+geom_histogram()+theme_bw()
ggplot(data=db, aes(x=caf))+geom_histogram()+theme_bw()
ggplot(data=dnb, aes(x=caf))+geom_histogram()+theme_bw()