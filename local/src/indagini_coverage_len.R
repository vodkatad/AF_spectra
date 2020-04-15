library(gridExtra)
library(ggplot2)

setwd("/mnt/trcanmed/snaketree/prj/AF_spectra/dataset/CRC1307_platypus_nobin")
d <- read.table("CRC1307LMO-0-B.callable.bed.gz",sep="\t", header=F)
colnames(d) <- c("chr", "b","e", "cn")
d$len <- d$e-d$b


#filelen <- "../../local/share/data/CRC1307_clones_mutect/wgs_calling_regions.hg38.bed.gz"
filelen <- "../../local/share/data/CRC1307_clones_mutect/callable_covered.bed.gz"
lengths <- read.table(gzfile(filelen), sep="\t", header=F, stringsAsFactors = F)
#other callable_covered.bed.gz TODO proporzioni
colnames(lengths) <- c("chr", "b","e", "x","y","z")
lengths$len <- lengths$e-lengths$b
lengths <- lengths[!lengths$chr %in% c("chrY","chrX"),]

#> sum(lengths$len)
#[1] 2745214431
#yabbadabbadou!

singlechr <- function(chr, data, len) {
  bedcn <- data[data$chr == chr,]
  noxnoy <- sum(len[len$chr==chr, 'len'])
  #noxnoy <- 2745214431 #noy <- 29001108840  was it on all? probably yes, but it seems reasonable to work only on the callable regions
  l <- unique(data$cn)
  allcn <- l[order(l)]
  allcn <- allcn[seq(1,6)]
  sums <- sapply(allcn, function(x) sum(bedcn[bedcn$cn <= x,"len"]))
  data <- data.frame(perc=sums/noxnoy,cn=allcn)
  return(ggplot(data = data, aes(x = cn, y = perc)) + geom_line() + geom_point()+theme_bw()+ggtitle(paste0("cum genome ", chr)))
}

plots <- lapply(unique(d$chr), singlechr, d, lengths)

#solo il 4 ha 8% di cn==1

#pdf("coverage_vs_calling_regions.pdf"); 
grid.arrange(grobs=plots[c(4,seq(8,13))]);
#dev.off()


# comparing covered vs callable
filelen1 <- "../../local/share/data/CRC1307_clones_mutect/wgs_calling_regions.hg38.bed.gz"
alengths <- read.table(gzfile(filelen1), sep="\t", header=F, stringsAsFactors = F)
#other callable_covered.bed.gz TODO proporzioni
colnames(alengths) <- c("chr", "b","e", "x","y","z")
alengths$len <- alengths$e-alengths$b
alengths <- alengths[!alengths$chr %in% c("chrY","chrX"),]

cov20x <- "../../local/share/data/CRC1307_clones_mutect/callable_20x.covered.bed.gz"
cov20 <- read.table(gzfile(cov20x), sep="\t", header=F, stringsAsFactors = F)
#other callable_covered.bed.gz TODO proporzioni
colnames(cov20) <- c("chr", "b","e", "x","y","z")
cov20$len <- cov20$e-cov20$b
cov20 <- cov20[!cov20$chr %in% c("chrY","chrX"),]

chrs <- unique(alengths$chr)
percs <- sapply(chrs, function(x) { 
  callable <- sum(alengths[alengths$chr==x,"len"])
  covered <- sum(lengths[lengths$chr==x,"len"])
  c20 <- sum(cov20[cov20$chr==x,"len"])
  
  #covered/callable
  return(c(callable, covered, c20))
  })


d <- t(percs)
colnames(d) <- c("callable","covered1xQ20","covered20xQ0")
library(reshape2)
library(scales)
long <- melt(d)
colnames(long) <- c("chr","type","len")
ggplot(long,aes(x = chr, y = len, fill = type)) +geom_bar(position = "dodge",stat = "identity")+theme_bw()+theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("Coverage 1x in all CRC1307 samples vs callable genome")
