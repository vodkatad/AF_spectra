#!/usr/bin/env Rscript
library(corrplot)
library(ggplot2)

subcl_f <- snakemake@input[['subcl']]
mr_f <- snakemake@params[['mr']]
outplot_f <- snakemake@output[['plot']]


subcl <- read.table(subcl_f, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(subcl) <- c('sample', 'n', 'higherAF', 'lowerAF')
mr <- read.table(mr_f, sep="\t", stringsAsFactors=FALSE, header=TRUE)

subcl$end <- sapply(subcl$sample, function(x) {y<-strsplit(x, '.', fixed=TRUE)[[1]][1]; return(y[1])})

save.image('pi.Rdata')
# remove in vivo
mr <- mr[!mr$class %in% c('MA', 'MC', 'MI'),]

all_thr <- unique(subcl[, c('lowerAF', 'higherAF')])
lowers <- sort(unique(all_thr$lowerAF))
highers <- sort(unique(all_thr$higherAF))
cors <- data.frame(matrix(NA, nrow=length(lowers), ncol=length(highers)))
rownames(cors) <- round(lowers, digits=2)
colnames(cors) <- round(highers, digits=2)
for (l_i in seq(1, length(lowers))) {
  for (h_i in seq(1, length(highers))) {
    if (lowers[l_i] < highers[h_i]) {
      data <- subcl[subcl$lowerAF == lowers[l_i] & subcl$higherAF == highers[h_i],]
      m <- merge(data, mr, by="end")
      cors[l_i, h_i] <- cor(m$n, m$MR_EDU)
    }
  }
}


png(outplot_f)
corrplot(as.matrix(cors))
graphics.off()

q(status=0, save="no")

l_i <- 2
h_i <- 5
lowers[l_i]
highers[h_i]
data <- subcl[subcl$lowerAF == lowers[l_i] & subcl$higherAF == highers[h_i],]
m <- merge(data, mr, by="end")
m$nn <- m$n / m$len_cnok
ggplot(data=m, aes(x=log10(MR_EDU), y=nn))+geom_point()+theme_bw()


l_i <- 1
h_i <- 1
lowers[l_i]
highers[h_i]
data <- subcl[subcl$lowerAF == lowers[l_i] & subcl$higherAF == highers[h_i],]
m <- merge(data, mr, by="end")

ggplot(data=m, aes(x=log10(MR_EDU), y=n))+geom_point()+theme_bw()
