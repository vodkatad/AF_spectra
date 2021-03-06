#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
mutfile <- args[1]
hetmuts <- args[2]
hascn <- as.logical(args[3])
output <- args[4]
pthr <- as.numeric(args[5])
log <- args[6]

lookfor <- read.table(gzfile(hetmuts), header=FALSE, stringsAsFactors=FALSE)
colnames(lookfor) <- c("id","homet","refreads","mutreads","af", "where")
muts <- read.table(gzfile(mutfile), header=FALSE, stringsAsFactors=FALSE)

total_gs_het <- nrow(lookfor)
# now we look for the muts in muts that are compatible with an AF of 0.5 with a binomial
# test, correcting for cn if it's there.
#(base) data@clo:~/work/AF_spectra/dataset/mixology$ zcat CRC0277LMX0A03202TUMD02000.nofilter.tsv.gz 
#chr18:22982188:G:A	1/1	0	220	0.996
if (!hascn) {
  colnames(muts) <- c("id","homet","refreads","mutreads","af")
  muts$cn <- 1
} else {
  colnames(muts) <- c("id","homet","refreads","mutreads","af","cn")
}

rbinom <- function(mut) {
  refreads <- mut[1]
  mutreads <- mut[2]
  cn <- mut[4]
  if (mutreads == 0) { # a bunch of calls for biod..
    return(FALSE)
  }
  tot <- refreads+mutreads
  pvals <- c()
  #pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
  if (cn == 1) {
    p <- 0.25
    bin <- binom.test(mutreads, tot, p=p, alternative="greater")
    pvals <- bin$p.value
  } else {
    # if cn == 3 we can be "het" being 1/3 or 2/3, if cn == 4 ...
    stop("CN != 1 still not implemented")
  }
  return(any(pvals < pthr))
}

consideredhet <- apply(muts[, c(3,4,5,6)], 1, rbinom)

save.image("pippo.RData")
hets <- muts[consideredhet,]
total_found_het <- nrow(hets)

common <- intersect(hets$id, lookfor$id)
lost <- setdiff(lookfor$id, hets$id)
lost_nocalled <- setdiff(lookfor$id, muts$id)
callable <- intersect(lookfor$id, muts$id)
lost_binomial <- setdiff(callable, hets$id)
fromnowhere <- setdiff(hets$id, lookfor$id)
callsfromnowhere <- setdiff(muts$id, lookfor$id)

info <- data.frame(what=c("total_gs_het","het_found","common","not_found", "not_found_call", "not_found_binomial", "appeared_binomial", "appeared_call"), n=c(total_gs_het, total_found_het, length(common), length(lost), length(lost_nocalled), length(lost_binomial), length(fromnowhere), length(callsfromnowhere)))

muts$het <- "no"
muts[consideredhet,]$het <- "yes"

muts$gs_het <- "no"
muts[muts$id %in% lookfor$id,]$gs_het <- "yes"

write.table(muts, gzfile(output), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(info, log, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
