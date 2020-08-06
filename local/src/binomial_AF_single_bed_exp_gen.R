#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
mutfile <- args[1]
hascn <- as.logical(args[2])
cns <- args[3]
pthr <- as.numeric(args[4])
output <- args[5]
logfile <- args[6]
outputdir <- args[7]

if (cns != "2,3") {
    stop("accepted_cn right now are hard-coded for 2,3")
}

accepted_cn <- as.numeric(unlist(strsplit(cns, ',', fixed=TRUE)))

muts <- read.table(mutfile, header=FALSE, stringsAsFactors=FALSE)
#chr1    181112  181113  chr1:181113:A:G:12:25:3
colnames(muts) <- c("chr","b","e","chrbis","b_1based", "ref", "alt", "refreads", "altreads", "af", "cn")
#R is an eternal source of misteries and painful data transformations
#https://stackoverflow.com/questions/15618527/why-does-as-matrix-add-extra-spaces-when-converting-numeric-to-character
rownames(muts) <- apply(sapply(muts, format, trim = TRUE, justify="none"), 1, function(x) {paste(x[1], x[5], x[6], x[7], sep=":")})

rbinom <- function(mut, generation) {
  refreads <- mut[1]
  mutreads <- mut[2]
  cn <- mut[3]
  if (mutreads == 0) { # a bunch of calls for biod..
    cat('caacca')
    exit(1)
  }
  tot <- refreads+mutreads
  if (cn == 2) { 
    p1 <- limits2[[generation]][1]
    p2 <- limits2[[generation]][2]
    bin1 <- binom.test(mutreads, tot, p=p1, alternative="less")
    bin2 <- binom.test(mutreads, tot, p=p2, alternative="greater")
    pval1 <- bin1$p.value
    pval2 <- bin2$p.value
  } else if (cn == 3) {
    p1 <- limits3[[generation]][1]
    p2 <- limits3[[generation]][2]
    bin1 <- binom.test(mutreads, tot, p=p1, alternative="less")
    bin2 <- binom.test(mutreads, tot, p=p2, alternative="greater")
    pval1 <- bin1$p.value
    pval2 <- bin2$p.value
  }
  return(c(pval1, pval2))
}

save.image("pippo.RData")

# limits up to 4 generations for cn2 and cn3
# 1/4 1/8 1/16 1/32
# lower is always 1/2 / 1/3 as we want to be inclusive
# "lower" upper limit for CN=3 is ok only as 1/3 (not 2/3 too):
# bino <- function(mutreads, tot) {
#   p1 <- binom.test(mutreads, tot, p=1/3, alternative="less")$p.value
#   p2 <- binom.test(mutreads, tot, p=2/3, alternative="less")$p.value
#   return((list(one=p1<0.05, two=p2<0.05)))
# }
# mut <- seq(1, 100)

# res <- lapply(mut, bino, 100)
# check <- res[sapply(res, function(x) isTRUE(x$one) & isFALSE(x$two)) ]
limits2 <- list(c(1/2,1/8), c(1/2,1/16), c(1/2,1/32), c(1/2, 1/64))
limits3 <- list(c(1/3,1/12), c(1/3,1/24), c(1/3,1/48), c(1/3,1/96))
# 1/2 - 1/128
# 1/3 - 1/192
muts$targetcn <- ifelse(muts$cn %in% accepted_cn, 1, 0)
muts <- muts[muts$targetcn==1,]
res <- NULL
log <- data.frame(generation=seq(1,4))
for (i in 1:4) {
  bires <-  t(apply(muts[, c(8,9,11)], 1, rbinom, i))
  muts[,"binomp1"] <- bires[,1]
  muts[,"binomp2"] <- bires[,2]
  muts$pBH1 <- p.adjust(muts$binomp1, method="BH")
  muts$pBonf1 <- p.adjust(muts$binomp1, method="bonferroni")
  muts$pBH2 <- p.adjust(muts$binomp2, method="BH")
  muts$pBonf2 <- p.adjust(muts$binomp2, method="bonferroni")
  muts$generation <- i
  if (is.null(res)) {
    res <- muts
  } else {
    res <- rbind(res, muts)
  }
  log[log$generation == i,"total"] <- nrow(muts)
  log[log$generation == i,"sign"] <- nrow(muts[muts$binomp1 < pthr & muts$binomp2 < pthr,])
  log[log$generation == i,"signbonf"] <- nrow(muts[muts$pBonf1 < pthr & muts$pBonf2 < pthr,])
  log[log$generation == i,"signbh"] <- nrow(muts[muts$pBH1 < pthr & muts$pBH2 < pthr,])
}

histo <- function(d, column, cn) {
    column1 <- paste0(column, '1')
    column2 <- paste0(column, '2')
    d <- d[d$cn == cn,]
    tot <- nrow(d)
    sign <- d[d[, column1] < pthr & d[,column2] < pthr,]
    ggplot(sign, aes(x=altreads)) + geom_histogram(binwidth=1)+xlab('N. mutated reads')+ylab('N.muts')+theme_bw()+ggtitle(paste('Total m: ', tot, 'Called m:', nrow(sign)))
    ggsave(paste0(outputdir,'/histo_', column, '_', cn,'.png'))
}   

allcn <- c(2,3)
allp <- c('binomp','pBH','pBonf')
wanted <- expand.grid(allcn,allp)

garbage <- apply(wanted, 1, function(x) { histo(res, x[2], x[1]) })


write.table(res, gzfile(output), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(log, logfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
