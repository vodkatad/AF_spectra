mutfile <- '/scratch/trcanmed/AF_spectra/dataset/CRC1307_mutect_nobin/CRC1307-09-0.var_cnv.tsv.gz'
muts <- read.table(gzfile(mutfile), header=FALSE, stringsAsFactors=FALSE)
colnames(muts) <- c("chr","b","e","id")
muts$cn <- sapply(muts$id, function(x) { strsplit(x, ":")[[1]][8]})
muts$cn <- as.numeric(muts$cn)
maxcn <- 4
muts <- muts[muts$cn <= maxcn,]

muts$af <- as.numeric(sapply(muts$id, function(x) { strsplit(x, ":")[[1]][7]}))
ggplot(data=muts, aes(x=af)) + geom_histogram(size=0.8, bins=60)+theme_bw()+theme(text=element_text(size=20))+facet_wrap(~cn)

o <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1307/platypus_nobin/o', sep="\t")

oo <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1307/platypus_nobin/oo', sep="\t")

distances <- function(chr, data) {
  dd <- data[data$V2=="chr1",]$V3
  af <- data[data$V2=="chr1",]$V11
  dd2 <- dd[-1]
  dd <- dd[-length(dd)]
  data.frame(dist=dd2-dd, af=af[-1])
}

meh <- lapply(levels(o$V2), distances, o )
dists <- do.call(rbind, meh)


cmeh <- lapply(levels(oo$V2), distances, oo )
cdists <- do.call(rbind, cmeh)

pdata <- data.frame(dist=c(cdists$dist, dists$dist), clone=c(rep('c', nrow(cdists)), rep('e', nrow(dists))), af=c(cdists$af, dists$af))
pdata$ldist <- log2(pdata$dist)
pdata$ldist <- cut(pdata$ldist ,breaks=c(0,5,20, max(pdata$ldist)))
ggplot(pdata, aes(x=log2(dist), color=clone))+geom_density()

ggplot(pdata, aes(x=af, fill=ldist))+geom_histogram(position="dodge")+theme_bw()+facet_wrap(~clone)
