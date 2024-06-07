#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
bedcnfile <- args[1]
cns <- args[2]
outputbarplot <- args[3]
outputcumplot <- args[4]
outputlen <- args[5]
outputlen2 <- args[7]

accepted_cn <- as.numeric(unlist(strsplit(cns, ',', fixed=TRUE)))

bedcn <- read.table(gzfile(bedcnfile), header=FALSE, stringsAsFactors=FALSE)
colnames(bedcn) <- c("chr","b","e","cnmy","cnthey")
bedcn$len <- bedcn$e - bedcn$b

#wantedcn_callable <- bedcn[bedcn$cn %in% accepted_cn,]
#wantedcn_callable <- bedcn[bedcn$cnmy %in% accepted_cn & bedcn$cnthey %in% accepted_cn,]
wantedcn_callable <- bedcn

write.table(wantedcn_callable, gzfile(args[6]), sep="\t", quote=F, row.names=F, col.names=F)

len <- sum(wantedcn_callable$len * wantedcn_callable$cnmy)
sink(outputlen)
cat(len)
cat("\n")
sink()

l <- unique(bedcn$cnmy)
lens <- sapply(l, function(x) { sum(bedcn[bedcn$cnmy==x,"len"])})
#pl <- data.frame(len=lens, cn=l)
nfrags <- sapply(l, function(x) { length(bedcn[bedcn$cnmy==x,"len"])})
#mean <-sapply(l, function(x) { mean(d[d$cn==x,"len"])})
pl <- data.frame(totlen=lens, cn=l, n=nfrags)
ggplot(data=pl, aes(x=cn, y=totlen, fill=nfrags)) +geom_bar(stat="identity")+scale_fill_gradient2(low="blue", mid="white", high="red", space="Lab")
ggsave(outputbarplot)

#(align_recalibrate) [CRC1307]egrassi@hactarlogin$ zcat /home/egrassi/bit/task/annotations/dataset/gnomad/wgs_calling_regions.hg38.bed.gz | filter_1col 1 /home/egrassi/WGS/local/share/data/CRC1307_simul/chrs | bawk '{print $3-$2}' > o
#> o <- read.table("o")
#[1] 2745214431
#
noxnoy <- 2745214431 #noy <- 29001108840  was it on all? probably yes, but it seems reasonable to work only on the callable regions
allcn <- l[order(l)]
sums <- sapply(allcn, function(x) sum(bedcn[bedcn$cnmy <= x,"len"]))
data <- data.frame(perc=sums/noxnoy,cn=allcn)
ggplot(data = data, aes(x = cn, y = perc)) + geom_line() + geom_point()+theme_bw()+ggtitle("cumulative genome covered")
ggsave(outputcumplot)

save.image('p.Rdata')
sink(outputlen2, append=T)
write.table(data, sep="\t", quote=FALSE, row.names=FALSE)
cat("\n")
sink()
