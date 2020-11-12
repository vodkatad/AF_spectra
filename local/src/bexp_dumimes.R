d <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1502_150x/all_bexp_different_coverage.tsv', sep="\t", header=TRUE)
d$coverage <- factor(d$coverage, levels=paste0(c(seq(10,100, by=10), 150, 200), 'x'))
n_divisions <- list('1'=1, '2'=3, '3'=7,'4'=15)

d$rate <- rep(0, nrow(d))
for (i in seq(1, nrow(d))) {
  d[i,]$rate <- d[i,]$called_in_gen / (n_divisions[[d[i,]$generation]] * d[i,]$len)
}
estimate <- 1.42820176969739e-09

d1 <- d[d$sample=="CRC1502-03-0",]
ggplot(data=d1, aes(x=coverage, y=log10(rate), shape=correction))+geom_point()+theme_bw()+facet_wrap(~generation)+geom_hline(yintercept=log10(estimate), color='red')+ggtitle(unique(d1$sample))

d1 <- d[d$sample=="CRC1502-03-1-A",]
ggplot(data=d1, aes(x=coverage, y=log10(rate), shape=correction))+geom_point()+theme_bw()+facet_wrap(~generation)+geom_hline(yintercept=log10(estimate), color='red')+ggtitle(unique(d1$sample))


d <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1307_150x/all_bexp_different_coverage.tsv', sep="\t", header=TRUE)
d$coverage <- factor(d$coverage, levels=paste0(c(seq(10,100, by=10), 150, 200), 'x'))
n_divisions <- list('1'=1, '2'=3, '3'=7,'4'=15)

d$rate <- rep(0, nrow(d))
for (i in seq(1, nrow(d))) {
  d[i,]$rate <- d[i,]$called_in_gen / (n_divisions[[d[i,]$generation]] * d[i,]$len)
}

estimate <- 2.56764427616427e-09
d1 <- d[d$sample=="CRC1307-02-1-E",]
ggplot(data=d1, aes(x=coverage, y=log10(rate), shape=correction))+geom_point()+theme_bw()+facet_wrap(~generation)+geom_hline(yintercept=log10(estimate), color='red')+ggtitle(unique(d1$sample))

d1 <- d[d$sample=="CRC1307-09-1-B",]
estimate <- 2.39506626830393e-09
ggplot(data=d1, aes(x=coverage, y=log10(rate), shape=correction))+geom_point()+theme_bw()+facet_wrap(~generation)+geom_hline(yintercept=log10(estimate), color='red')+ggtitle(unique(d1$sample))
