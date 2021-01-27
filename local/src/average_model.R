#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]

our <- read.table(infile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(our) <- c('sample','estimate')
our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})

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
ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'estimate'], LEVEL) })
colnames(ic_clones) <- unique(our$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)
write.table(pdata, file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)