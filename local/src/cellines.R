library(ggplot2)
d <- read.table('/tmp/o')
colnames(d) <- c('filename', 'npass')
d$sample <- gsub('.+_BASELINE__', '', d$filename)
d$sample <- gsub('_BULK.pass.vcf.gz', '', d$sample, fixed=T)


d<- d[order(d$sample),]
d$line <- ifelse(grepl('DIFI', d$sample), 'DIFI', 'WIDR')
ggplot(data=d, aes(x=sample, y=npass, fill=line))+geom_col()+theme_bw(base_size=15)+scale_y_sqrt(n.breaks=11)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#GRCh38.d1.vd1__
d <- read.table('/tmp/oo')
colnames(d) <- c('filename', 'npass')
d$sample <- gsub('GRCh38.d1.vd1__', '', d$filename, fixed=T)
d$sample <- gsub('_B.+.pass.vcf.gz', '', d$sample, fixed=F)


d<- d[order(d$sample),]
d$line <- ifelse(grepl('DIFI', d$sample), 'DIFI', 'WIDR')
ggplot(data=d, aes(x=sample, y=npass, fill=line))+geom_col()+theme_bw(base_size=15)+scale_y_sqrt(n.breaks=11)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


dd <- read.table('/mnt/trcanmed/snaketree/data/MA_cell_lines_MR/to_MA/depmap_woverlap.tsv', sep="\t", header=T)
dd$sample <- gsub('_BULK', '', dd$sample, fixed=T)
dd$sample <- gsub('_BASELINE', '', dd$sample, fixed=F)
dd$frac <- dd$n_overlap/dd$n_depmap
dd$line <- ifelse(grepl('DIFI', dd$sample), 'DIFI', 'WIDR')

ggplot(data=dd, aes(x=sample, y=frac, fill=line))+geom_col()+theme_bw(base_size=15)+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
