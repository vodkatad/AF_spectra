#egrassi@godot:/scratch/trcanmed/AF_spectra/datasetV2$ zcat */mutect_nobin/*.ovcnokdelta.tsv.gz > o
#egrassi@godot:/scratch/trcanmed/AF_spectra/datasetV2$ zcat CRC1307*/platypus_nobin_00/*.ovcnokdelta.tsv.gz CRC1502*/platypus_nobin_00/*.ovcnokdelta.tsv.gz  CRC0441*/platypus_nobin_00/*.ovcnokdelta.tsv.gz > oo
library(ggsignif)
library(ggplot2)

d <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/o', sep="\t", header=F)
colnames(d) <- c('id', 'chr', 'b', 'e', 'chr2', 'b2', 'ref', 'alt', 'tot', 'mut', 'AF', 'cn', 'n1', 'n2', 'n3', 'class')
ggplot(data=d, aes(x=class, y=AF))+geom_violin()+geom_boxplot(width=0.1)+theme_bw()+ 
  geom_signif(comparisons = list(c("gain", "loss")), map_signif_level = FALSE)


dd <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/oo', sep="\t", header=F)
colnames(dd) <- c('id', 'chr', 'b', 'e', 'chr2', 'b2', 'ref', 'alt', 'tot', 'mut', 'AF', 'cn', 'n1', 'n2', 'n3', 'class')
ggplot(data=dd, aes(x=class, y=AF))+geom_violin()+geom_boxplot(width=0.03)+theme_bw()+ 
  geom_signif(comparisons = list(c("gain", "loss")), map_signif_level = FALSE)


d <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/mutect_nobin/o', sep="\t", header=F)
colnames(d) <- c('id', 'chr', 'b', 'e', 'chr2', 'b2', 'ref', 'alt', 'tot', 'mut', 'AF', 'cn', 'n1', 'n2', 'n3', 'class')
ggplot(data=d, aes(fill=class, x=mut))+geom_bar(position='dodge')+theme_bw()
table(d[d$mut==1,'class'])
table(d[d$mut==2,'class'])

d <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/platypus_nobin_00/o', sep="\t", header=F)
colnames(d) <- c('id', 'chr', 'b', 'e', 'chr2', 'b2', 'ref', 'alt', 'tot', 'mut', 'AF', 'cn', 'n1', 'n2', 'n3', 'class')
ggplot(data=d, aes(fill=class, x=mut))+geom_bar(position='dodge')+theme_bw()
table(d[d$mut==1,'class'])
table(d[d$mut==2,'class'])
