library(ggplot2)

n_f  <- snakemake@input[['n']]
genome_f  <- snakemake@input[['gb']]
ave_MR_f  <- snakemake@input[['ave_MR']]
#colors <- snakemake@input[['palette']]

#log_f <- snakemake@log[['log']]
estim_age_f <- snakemake@output[['estim_age']]

save.image(paste0(estim_age_f, '.Rdata'))

d <- read.table(n_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(d) <- c('#_SNVs', 'sample')

d$type <- ifelse(grepl('bulk', d$sample), 'parental', 'T0_clone')
d$model <- rep('', nrow(d))
d[d$type=="parental", 'model'] <- sapply(strsplit(d[d$type=="parental", 'sample'], "_"), '[[', 2)
d[d$type!="parental", 'model'] <- sapply(strsplit(d[d$type!="parental", 'sample'], "-"), '[[', 1)
#d$model <- ifelse(d$type=="parental", sapply(strsplit(d$sample, "_"), '[[', 2), sapply(strsplit(d$sample, "-"), '[[', 1))
d$model <-  paste0(d$model, ifelse(!grepl('\\d$', d$model), '', ifelse(d$model=="CRC0282", 'PR', 'LM')))

gen_len <- read.table(genome_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(gen_len) <- c('sample', 'Gb')

# average len for each model T1. This is the len of T1-T0.ovlength, on which we computed the MR
gen_len <- gen_len[grepl('-1-', gen_len$sample),]
gen_len$model <- sapply(gen_len$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
gen_len$clone <- sapply(gen_len$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
gen_len$clone2 <- sapply(gen_len$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
gen_len$model_clone <- paste0(gen_len$model, "_", gen_len$clone)

gen_len$model <- paste0(gen_len$model, ifelse(!grepl('\\d$', gen_len$model), '', ifelse(gen_len$model=="CRC0282", 'PR', 'LM')))

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
ic_clones <- sapply(unique(gen_len$model), function(x) { confidence_interval(gen_len[gen_len$model==x,'Gb'], LEVEL) })
colnames(ic_clones) <- unique(gen_len$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)


ave_MR <- read.table(ave_MR_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
ave_MR$model <- rownames(ave_MR)

ave_gb <- pdata
# ave t0 n muts
n_clones <- d[d$type== "T0_clone",]
LEVEL <- 0.99
ic_clones <- sapply(unique(n_clones$model), function(x) { confidence_interval(n_clones[n_clones$model==x,'#_SNVs'], LEVEL) })
colnames(ic_clones) <- unique(n_clones$model)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

ave_n <- pdata

bulk <- d[d$type== "parental",]
# Now we have average T0 n muts in ave_n, average genome len in ave_gb, ave MR in ave_MR and bulk/parental Ns in bulk
colnames(ave_n)[3] <- 'avg_n'
colnames(ave_gb)[3] <- 'avg_gb'
colnames(ave_MR)[3] <- 'avg_MR'

m <- merge(ave_n, ave_gb, by="model")
m2 <- merge(m, ave_MR, by="model")
m3 <- merge(m2, bulk, by="model")

m3[,grepl('lower', colnames(m3))] <- NULL
m3[,grepl('upper', colnames(m3))] <- NULL
m3$xmodel <- NULL
m3$sample <- NULL
m3$type <- NULL
m3$avg_gb <- m3$avg_gb / 10**9
m3$age <- (m3$avg_n  - m3$`#_SNVs`)/(m3$avg_gb * m3$avg_MR)

write.table(m3, file=estim_age_f, row.names=FALSE, quote=FALSE, sep="\t")
