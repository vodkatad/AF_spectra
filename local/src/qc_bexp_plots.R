d <- read.table("/scratch/trcanmed/AF_spectra/dataset/CRC1307_150x/mutect_nobin_100x/CRC1307-02-1-E.bexpcalls.gen.tsv.gz", sep="\t", header=T)
head(d)
library(ggplot2)
dcn2 <- d[d$cn==2 & d$generation==4,]
#af <- dcn2$af; excum <- sapply(1:length(af),function(i) sum(af[i] <= af[1:length(af)]))
dim(dcn2)
af <- dcn2$af
excum <- sapply(1:length(af),function(i) sum(af[i] <= af[1:length(af)]))
data <- data.frame(af=af, cum=excum)
dim(data)
#+ coord_trans(y = "log10")
ggplot(data, aes(x=log10(af), y=log10(excum))) + geom_point()+geom_smooth(method=lm)+theme_bw()+xlab('log10(AF)')+ylab('Log10 Cumulative n. of muts > AF')
head(dcn2)
dsign <- dcn2[dcn2$pBH2 < 0.05 & dcn2$pBH1 < 0.05,]
af <- dsign$af
excum <- sapply(1:length(af),function(i) sum(af[i] <= af[1:length(af)]))
data <- data.frame(af=af, cum=excum)
dim(data)

ggplot(data, aes(x=log10(af), y=log10(excum))) + geom_point()+theme_bw()+xlab('log10(AF)')+ylab('Log10 Cumulative n. of muts > AF')+geom_smooth(method=lm)

# sottoriva
ddd <- dsign
exsubcl <- ddd$af
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
invf <- 1/exsubcl
maxi <- length(invf)
labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
model <- lm(excum~invf)
sfit <- summary(model)
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit$r.squared, digits=3)))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model, col="red")


## bonf
dsign <- dcn2[dcn2$pBonf1 < 0.05 & dcn2$pBonf2 < 0.05,]
af <- dsign$af
excum <- sapply(1:length(af),function(i) sum(af[i] <= af[1:length(af)]))
data <- data.frame(af=af, cum=excum)

ggplot(data, aes(x=log10(af), y=log10(excum))) + geom_point()+theme_bw()+xlab('log10(AF)')+ylab('Log10 Cumulative n. of muts > AF')+geom_smooth(method=lm)

# sottoriva
ddd <- dsign
exsubcl <- ddd$af
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
invf <- 1/exsubcl
maxi <- length(invf)
labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
model <- lm(excum~invf)
sfit <- summary(model)
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit$r.squared, digits=3)))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model, col="red")


sottoriva <- function(data) {
  ddd <- data
  exsubcl <- ddd$af
  excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
  invf <- 1/exsubcl
  maxi <- length(invf)
  labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
  model <- lm(excum~invf)
  sfit <- summary(model)
  plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit$r.squared, digits=3)))
  oi <- invf[order(invf)]
  oex <- exsubcl[order(-exsubcl)]
  axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
  abline(model, col="red")
  return(sfit)
}

#### 
path30x <- '/scratch/trcanmed/AF_spectra/dataset/CRC1502/mutect_nobin/CRC1502-03-1-A.var_cnv.tsv.gz'
path120x <- '/scratch/trcanmed/AF_spectra/dataset/CRC1502_150x/mutect_nobin_30x/CRC1502-03-1-A.var_cnv.tsv.gz'

load <- function(path) {
  d <- read.table(path, sep="\t", stringsAsFactors=FALSE, header=FALSE)
  exploded <- lapply(d$V4, function(x) {strsplit(x, ":", fixed=TRUE)})
  colnames(d) <- c('chr', 'b','e','id')
  rownames(d) <- d$id
  d$id <- NULL
  d$af <- as.numeric(sapply(exploded, function(x) { x[[1]][7]}))
  d$cn <- as.numeric(sapply(exploded, function(x) { x[[1]][8]}))
  return(d)
}

lc <- load(path30x)
hc <- load(path120x)

s_lc <- sottoriva(lc[lc$cn==2 & lc$af > 0.12 & lc$af < 0.24,])
s_hc <- sottoriva(hc[hc$cn==2 & hc$af > 0.12 & hc$af < 0.24,])

s_lc_2 <- sottoriva(lc[lc$cn==2 & lc$af > 0.12 & lc$af < 0.4,]) # less than expected...
s_hc_2 <- sottoriva(hc[hc$cn==2 & hc$af > 0.12 & hc$af < 0.4,])

# need to fix the intercept?
#linear model with slope μ/β and intercept –μ/(β fmax). We exploited the intercept constraint to perform a more restrictive fit using the model y=m(x-1/fmax)+0.

bulk_path <- '/scratch/trcanmed/AF_spectra/dataset/second_shipment_bulk/mutect_nobin/CRC1502LMO-0-B.var_cnv.tsv.gz'
b <- load(bulk_path)
s_b <- sottoriva(b[b$cn==2 & b$af > 0.12 & b$af < 0.24,])

