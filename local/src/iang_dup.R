di <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset_IANG/CRCUECHPRO/platypus_nobin_00_allcn/CRCUECHPRO-13-1-D_CRCUECHPRO-13-0.ovcnokdelta.tsv.gz'), sep="\t", header=F)
ci <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset_IANG/CRCUECHPRO/platypus_nobin_00_allcn/CRCUECHPRO-13-1-C_CRCUECHPRO-13-0.ovcnokdelta.tsv.gz'), sep="\t", header=F)

di_g <- di[di$V16=="gain",]
ci_g <- ci[ci$V16=="gain",]

library(ggplot2)
ggplot(data=di_g, aes(x=V11))+geom_histogram()+theme_bw()
ggplot(data=ci_g, aes(x=V11))+geom_histogram()+theme_bw()

tdchr5 <- di_g[di_g$V2=="chr2",]
cchr5 <- ci_g[ci_g$V2=="chr2",]

ggplot(data=dchr5, aes(x=V11))+geom_histogram()+theme_bw()
ggplot(data=cchr5, aes(x=V11))+geom_histogram()+theme_bw()


dchr6 <- di_g[di_g$V2=="chr6",]
cchr6 <- ci_g[ci_g$V2=="chr6",]

ggplot(data=dchr6, aes(x=V11))+geom_histogram()+theme_bw()
ggplot(data=cchr6, aes(x=V11))+geom_histogram()+theme_bw()


dchr1 <- di_g[di_g$V2=="chr1",]
cchr1 <- ci_g[ci_g$V2=="chr1",]

ggplot(data=dchr1, aes(x=V11))+geom_histogram(aes(y=..ncount..))+theme_bw()+geom_vline(xintercept=1/4)
ggplot(data=cchr1, aes(x=V11))+geom_histogram(aes(y=..ncount..))+theme_bw()+geom_vline(xintercept=1/4)

dcn5 <- di_g[di_g$V12==5,]
ccn5 <- ci_g[ci_g$V12==5,]

ggplot(data=dcn5, aes(x=V11))+geom_histogram()+theme_bw()
ggplot(data=ccn5, aes(x=V11))+geom_histogram()+theme_bw()

ggplot(data=dcn5, aes(x=V11))+geom_histogram()+theme_bw()+geom_vline(xintercept=1/5)+
  geom_vline(xintercept=2/5)+geom_vline(xintercept=3/5)+geom_vline(xintercept=4/5)

ggplot(data=dcn5, aes(x=V11))+geom_histogram()+theme_bw()
ggplot(data=ccn5, aes(x=V11))+geom_histogram()+theme_bw()

ggplot(data=dcn5, aes(x=V11))+geom_histogram()+theme_bw()+geom_vline(xintercept=1/5)+
  geom_vline(xintercept=2/5)+geom_vline(xintercept=3/5)+geom_vline(xintercept=4/5)

dcn6 <- di_g[di_g$V12==6,]
ggplot(data=dcn6, aes(x=V11))+geom_histogram(bins=50)+theme_bw()+geom_vline(xintercept=1/6)+
  geom_vline(xintercept=2/6)+geom_vline(xintercept=3/6)+geom_vline(xintercept=4/6)+geom_vline(xintercept=5/6)

ccn3 <- ci_g[ci_g$V12==3,]
ggplot(data=ccn3, aes(x=V11))+geom_histogram()+theme_bw()+geom_vline(xintercept=1/3)+
  geom_vline(xintercept=2/3)+xlim(0,1)



dcn4 <- di_g[di_g$V12==4,]
ggplot(data=dcn4, aes(x=V11))+geom_histogram()+theme_bw()+geom_vline(xintercept=1/4)+
  geom_vline(xintercept=1/2)+geom_vline(xintercept=3/4)

dall <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset_IANG/CRCUECHPRO/platypus_nobin_00_allcn/CRCUECHPRO-13-1-D.calls.tsv.gz'), sep="\t", header=T, stringsAsFactors = F)

call <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset_IANG/CRCUECHPRO/platypus_nobin_00_allcn/CRCUECHPRO-13-1-C.calls.tsv.gz'), sep="\t", header=T, stringsAsFactors = F)

dall <- dall[dall$cn == 6,]
ggplot(data=dall, aes(x=af))+geom_histogram(bins=50)+theme_bw()+geom_vline(xintercept=2/3)+
  geom_vline(xintercept=1/3)+  geom_vline(xintercept=1/6)



call <- call[call$cn == 3,]
ggplot(data=call, aes(x=af))+geom_histogram()+theme_bw()+geom_vline(xintercept=2/3)+
  geom_vline(xintercept=1/3)

### 
dall$tot <- dall$refreads+dall$altreads
summary(dall$tot)
dall <- dall[dall$tot >= 37,]
ggplot(data=dall, aes(x=af))+geom_histogram(bins=50)+theme_bw()+geom_vline(xintercept=2/3)+
  geom_vline(xintercept=1/3)+  geom_vline(xintercept=1/6)

dcn6$tot <- dcn6$V9+dcn6$V10
summary(dcn6$tot)
dcn6 <- dcn6[dcn6$tot >= 34,]

ggplot(data=dcn6, aes(x=V11))+geom_histogram(bins=50)+theme_bw()+geom_vline(xintercept=1/6)+
  geom_vline(xintercept=2/6)+geom_vline(xintercept=3/6)+geom_vline(xintercept=4/6)+geom_vline(xintercept=5/6)


dall <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset_IANG/CRCUECHPRO/platypus_nobin_00_allcn/CRCUECHPRO-13-1-D.calls.tsv.gz'), sep="\t", header=T, stringsAsFactors = F)
dall2 <- dall[dall$chr=="chr2",]

ggplot(data=dall2, aes(x=af))+geom_histogram(bins=50)+theme_bw()+geom_vline(xintercept=1/5)
ggplot(data=dall2, aes(x=af))+geom_histogram(bins=50)+theme_bw()+geom_vline(xintercept=1/3)+geom_vline(xintercept=1/5)
