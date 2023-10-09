library(ggplot2)

n_f  <- snakemake@input[['n']]
ave_MR_f  <- snakemake@input[['ave_MR']]
#colors <- snakemake@input[['palette']]

#log_f <- snakemake@log[['log']]
estim_div_f <- snakemake@output[['estim_div']]

save.image(paste0(estim_div_f, '.Rdata'))

d <- read.table(n_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(d) <- c('sample', 'nSNVs', 'len' )

# average len for each model T1. This is the len of T1-T0.ovlength, on which we computed the MR
d$model <- sapply(d$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
d$clone <- sapply(d$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
d$clone2 <- sapply(d$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
d$model_clone <- paste0(d$model, "_", d$clone)

d$model <- paste0(d$model, ifelse(!grepl('\\d$', d$model), '', ifelse(d$model=="CRC0282", 'PR', 'LM')))

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- vec_sd
  
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}

LEVEL <- 0.99
ic_clones <- sapply(unique(d$model), function(x) { confidence_interval(d[d$model==x,'nSNVs'], LEVEL) })
colnames(ic_clones) <- unique(d$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

ave_n <- pdata

ic_clones <- sapply(unique(d$model), function(x) { confidence_interval(d[d$model==x,'len'], LEVEL) })
colnames(ic_clones) <- unique(d$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

ave_gb <- pdata

ave_MR <- read.table(ave_MR_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
ave_MR$model <- rownames(ave_MR)

# Now we have average T0 n muts in ave_n, average genome len in ave_gb, ave MR in ave_MR and bulk/parental Ns in bulk
colnames(ave_n)[3] <- 'avg_n'
colnames(ave_gb)[3] <- 'avg_gb'
colnames(ave_MR)[3] <- 'avg_MR'

m <- merge(ave_n, ave_gb, by="model")
m2 <- merge(m, ave_MR, by="model")

m2[,grepl('lower', colnames(m2))] <- NULL
m2[,grepl('upper', colnames(m2))] <- NULL
m2$xmodel <- NULL
m2$avg_gb <- m2$avg_gb / 10**9

m2$vivo_div <- m2$avg_n/(m2$avg_gb * m2$avg_MR)

write.table(m2, file=estim_div_f, row.names=FALSE, quote=FALSE, sep="\t")
