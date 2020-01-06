#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
mutfile <- args[1]
hascn <- as.logical(args[2])
cns <- args[3]
pthr <- as.numeric(args[4])
output <- args[5]

accepted_cn <- as.numeric(unlist(strsplit(cns, ',' fixed=TRUE)))

muts <- read.table(mutfile, header=FALSE, stringsAsFactors=FALSE)
#chr1    181112  181113  chr1:181113:A:G:12:25:3
colnames(muts) <- c("chr","b","e","chrbis","b_1based", "ref", "alt", "refreads", "altreads", "af", "cn")
rownames(muts) <- apply(muts, 1, function(x) {paste(x[1], x[5], x[6], x[7], sep=":")})

rbinom <- function(mut) {
  refreads <- mut[1]
  mutreads <- mut[2]
  cn <- mut[3]
  if (mutreads == 0) { # a bunch of calls for biod..
    return(FALSE)
  }
  tot <- refreads+mutreads
  pvals <- c()
  #pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
  if (cn == 1) {
    p <- 0.25
    bin <- binom.test(mutreads, tot, p=p, alternative="two.sided")
    pvals <- bin$p.value
  } else {
    # if cn == 3 we can be "het" being 1/3 or 2/3, if cn == 4 ...
    stop("CN != 1 still not implemented")
  }
  return(any(pvals < pthr))
}

consideredfounder <- apply(muts[, c(8,9,11)], 1, rbinom)

save.image("pippo.RData")
muts$founder <- 0
muts[consideredfounder,]$founder <- 1

muts$targetcn <- muts$cn %in% accepted_cn


write.table(muts, gzfile(output), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
