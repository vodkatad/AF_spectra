#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
mutfile <- args[1]
hascn <- as.logical(args[2])
cns <- args[3]
pthr <- as.numeric(args[4])
output <- args[5]
log <- args[6]

if (cns != "1,2,3") {
    stop("accepted_cn right now are hard-coded for 1,2,3")
}

accepted_cn <- as.numeric(unlist(strsplit(cns, ',', fixed=TRUE)))

muts <- read.table(mutfile, header=FALSE, stringsAsFactors=FALSE)
#chr1    181112  181113  chr1:181113:A:G:12:25:3
colnames(muts) <- c("chr","b","e","chrbis","b_1based", "ref", "alt", "refreads", "altreads", "af", "cn")
#R is an eternal source of misteries and painful data transformations
#https://stackoverflow.com/questions/15618527/why-does-as-matrix-add-extra-spaces-when-converting-numeric-to-character
rownames(muts) <- apply(sapply(muts, format, trim = TRUE, justify="none"), 1, function(x) {paste(x[1], x[5], x[6], x[7], sep=":")})

rbinom <- function(mut) {
  refreads <- mut[1]
  mutreads <- mut[2]
  cn <- mut[3]
  pval <- 1
  if (mutreads == 0) { # a bunch of calls for biod..
    return(pval)
  }
  tot <- refreads+mutreads
  #pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
  if (cn == 2) { # 1/(cn*2) could be the rule?
    p <- 0.5
    bin <- binom.test(mutreads, tot, p=p, alternative="less")
    pval <- bin$p.value
  } else if (cn == 3) {
    p <- 1/3
    bin <- binom.test(mutreads, tot, p=p, alternative="less")
    pval <- bin$p.value
  } else if (cn == 1) {
    p <- 1
    bin <- binom.test(mutreads, tot, p=p, alternative="less")
    pval <- bin$p.value
  }
  return(pval)
}

muts$targetcn <- ifelse(muts$cn %in% accepted_cn, 1, 0)
muts$binomp <- apply(muts[, c(8,9,11)], 1, rbinom)

save.image("pippo.RData")

muts$founder <- 0
muts$binomp <- p.adjust(muts$binomp, method="BH")
muts[muts$binomp < pthr,]$founder <- 1

#TODO LOG with some stats
nmut <- nrow(muts)
n_cnok <- sum(muts$targetcn)
n_cnko <- nrow(muts) - n_cnok
n_binomialok <- sum(muts$founder)
n_binomialko <- nrow(muts) - n_binomialok
x <- muts[muts$founder == 0 & muts$targetcn==1,]
n_cnok_binomialko <- nrow(x)
x <- muts[muts$founder == 1 & muts$targetcn==1,]
n_cnok_binomialok <- nrow(x)
x <- muts[muts$founder == 0 & muts$targetcn==0,]
n_cnko_binomialko <- nrow(x)
x <- muts[muts$founder == 1 & muts$targetcn==0,]
n_cnko_binomialok <- nrow(x)

info <- data.frame(what=c("nmut","n_cnok","n_cnko","n_binomialok","n_binomialko","n_cnok_binomialko","n_cnok_binomialok","n_cnko_binomialko","n_cnko_binomialok"), n=c(nmut,n_cnok,n_cnko,n_binomialok,n_binomialko,n_cnok_binomialko,n_cnok_binomialok,n_cnko_binomialko,n_cnko_binomialok))
write.table(info, log, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

write.table(muts, gzfile(output), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
