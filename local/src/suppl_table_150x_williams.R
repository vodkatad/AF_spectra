
williams_f  <- snakemake@input[['williams']]
MR_edu_f <- snakemake@input[['MR']]
bd_f  <- snakemake@input[['b_d']]
len_f <- snakemake@input[['gen_len']]
log_f <- snakemake@log[['log']]
outtsv <- snakemake@output[['out']]

our <- read.table(williams_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)

our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)

mredu <- read.table(MR_edu_f, sep="\t", header=FALSE)
colnames(mredu) <- c('sample', 'MRedu')
# Do we need to study williams only on CN 1-2-3?
# 1-b/d

bd <- read.table(bd_f, sep="\t", stringsAsFactors = FALSE, header=TRUE)
len <- read.table(len_f, sep="\t", stringsAsFactors = FALSE, header=TRUE)
colnames(bd)[1] <- 'sample'

m <- merge(our, len, by="sample")

bd$model <- sapply(bd$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
bd$clone <- sapply(bd$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
bd$model_clone <- paste0(bd$model, "_", bd$clone)

m2 <- merge(m, bd, by="model_clone")
m2$MRw <- (m2$intercept / m2$len)*(m2$X.b.d..b)

# MR_edu
#mredu$MRedu <- mredu$MRedu / 0.000000001
mm <- merge(mredu, m2, by="sample", by.y="sample.x")

res <- mm[, c('sample', 'r', 'intercept', 'len', 'b', 'b.d', 'MRw', 'MRedu')]

write.table(res, outtsv, quote=FALSE, sep="\t", row.names=FALSE)