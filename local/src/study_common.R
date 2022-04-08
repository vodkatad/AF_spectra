dA <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset/CRC1502/platypus_nobin_cn2/CRC1502-03-1-A.calls.tsv.gz'), sep="\t", header=TRUE)
dC <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset/CRC1502/platypus_nobin_cn2/CRC1502-03-1-C.calls.tsv.gz'), sep="\t", header=TRUE)
dD <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset/CRC1502/platypus_nobin_cn2/CRC1502-03-1-D.calls.tsv.gz'), sep="\t", header=TRUE)


n_found <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1502/platypus_nobin_cn2/CRC1502-03-0.sharedmuts.tsv', sep="\t", header=T)

m_dA <- merge(dD, n_found, by.y='ugains', by.x="row.names")

#> dim(dA)
#[1] 23916    14
#> dim(m_dA)
#[1] 694  16

#  CRC1502-03-1-A  CRC1502-03-0    1       691     1235003222      6.56589532532941e-09    1.05265413132127e-09    2.92938591792631e-09

# the not gained ones are removed for Cn filters in T0 sample?

g <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset/CRC1502/platypus_nobin_cn2/CRC1502-03-1-A_CRC1502-03-0.ovcnokdelta.tsv.gz'), sep="\t", header=F)

table(g[g$V1 %in% m_dA$Row.names,'V16'])

#common   gain   loss 
#0    691      0

#>  m_dA[! m_dA$Row.names %in% g$V1,]
#Row.names   chr        b        e chrbis b_1based ref alt refreads altreads        af cn targetcn binomp founder Freq
#250 chr17:21778905:A:G chr17 21778904 21778905  chr17 21778905   A   G       28        4 0.1428571  3        0      0       1    1
#280    chr18:77828:A:G chr18    77827    77828  chr18    77828   A   G        9        1 0.1111111  3        0      0       1    1
#368 chr22:16422510:A:G chr22 16422509 16422510  chr22 16422510   A   G       18        2 0.1111111  3        0      0       1    2

# ok!

m_dA <- m_dA[m_dA$cn == 2,]
table(m_dA$Freq)

ggplot(data=m_dA, aes(x=as.factor(Freq), y=af))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+theme_bw()
m_dA$cov <- m_dA$refreads + m_dA$altreads
m_dA$af2 <- m_dA$altreads / m_dA$cov
ggplot(data=m_dA, aes(x=as.factor(Freq), y=cov))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+theme_bw()
ggplot(data=m_dA, aes(x=as.factor(Freq), y=af2))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+theme_bw()
ggplot(data=m_dA, aes(x=as.factor(Freq), y=altreads))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+theme_bw()
wilcox.test(m_dA[m_dA$Freq==2, 'altreads'],m_dA[m_dA$Freq==1, 'altreads'])
wilcox.test(m_dA[m_dA$Freq==2, 'cov'],m_dA[m_dA$Freq==1, 'cov'])
# let's discuss quality
library('vcfR')
vcf <- read.vcfR('/scratch/trcanmed/AF_spectra/local/share/data/CRC1502/platypus_filtered.vcf.gz')

genoqual <- extract.gt(vcf, element = "GQ", as.numeric = TRUE) # ##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">

m_dA$id <- paste0(m_dA$chr, "_", m_dA$b+1)
mm <- m_dA[, c('Freq', 'id')]
m2 <- merge(mm, genoqual, by.x="id", by.y="row.names")

ggplot(data=m2, aes(x=as.factor(Freq), y=`CRC1502-03-1-D`))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+theme_bw()

# what are those:

#> m_dA[m_dA$altreads == 0,]
#Row.names   chr         b         e chrbis  b_1based ref alt refreads altreads af cn targetcn binomp founder Freq cov af2
#36   chr1:28354752:T:A  chr1  28354751  28354752   chr1  28354752   T   A       17        0  0  2        1      0       1    3  17   0
#138 chr15:20252465:G:A chr15  20252464  20252465  chr15  20252465   G   A       45        0  0  2        1      0       1    3  45   0
#197 chr16:54680859:A:G chr16  54680858  54680859  chr16  54680859   A   G       24        0  0  2        1      0       1    1  24   0
#234 chr17:16385615:C:T chr17  16385614  16385615  chr17  16385615   C   T       26        0  0  2        1      0       1    1  26   0
#244  chr17:7157284:T:C chr17   7157283   7157284  chr17   7157284   T   C       22        0  0  2        1      0       1    3  22   0
#245  chr17:7157288:C:G chr17   7157287   7157288  chr17   7157288   C   G       21        0  0  2        1      0       1    3  21   0
#509 chr4:167571484:T:A  chr4 167571483 167571484   chr4 167571484   T   A        7        0  0  2        1      0       1    3   7   0
#639 chr5:162488819:C:A  chr5 162488818 162488819   chr5 162488819   C   A       18        0  0  2        1      0       1    2  18   0


###INFO=<ID=MQ,Number=.,Type=Float,Description="Root mean square of mapping qualities of reads at the variant position">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant-quality/read-depth for this variant">

rms_mapq <- data.frame(row.names=row.names(genoqual),  mapq=extract.info(vcf, element = "MQ", as.numeric = TRUE)) # ##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">
#rms_mapq <- data.frame(row.names=row.names(genoqual),  mapq=extract.info(vcf, element = "QD", as.numeric = TRUE)) # ##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">
m2 <- merge(mm, rms_mapq, by.x="id", by.y="row.names")

ggplot(data=m2, aes(x=as.factor(Freq), y=mapq))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+theme_bw()

wilcox.test(m2[m2$Freq==2, 'mapq'],m2[m2$Freq==1, 'mapq'])

summary(m2[m2$Freq==2, 'mapq']); summary(m2[m2$Freq==1, 'mapq'])

# signal for QD:
#> summary(m2[m2$Freq==2, 'mapq']); summary(m2[m2$Freq==1, 'mapq'])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#10.01   20.82   22.95   22.43   24.50   44.05 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#10.28   17.22   20.34   20.70   23.77   49.10 
