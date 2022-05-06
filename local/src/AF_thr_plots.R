library(ggplot2)
data <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1502/AF_thr_MRedu', sep= "\t", header=F, stringsAsFactors = F)
colnames(data) <- c('file', 'clone', 'MR')
remove <- c('platypus_nobin_00/all.MR_ov', 'platypus_nobin_cn2/all.MR_ov', 'platypus_nobin_indels/all.MR_ov', 'platypus_nobin/all.MR_ov')
data <- data[!data$file %in% remove,]

data$thr <- gsub('platypus_nobin_',  '', data$file)
data$thr <- gsub('/all.MR_ov',  '', data$thr)
data[data$thr=="cn3", 'thr'] <- 0
data[data$thr=="003", 'thr'] <- 0.03
data[data$thr=="006", 'thr'] <- 0.06
data[data$thr=="01", 'thr'] <- 0.1
data[data$thr=="03", 'thr'] <- 0.3
data[data$thr=="05", 'thr'] <- 0.5
data[data$thr=="07", 'thr'] <- 0.7
data[data$thr=="09", 'thr'] <- 0.9


ggplot(data=data, aes(y=MR, color=clone, x=thr))+geom_point()+geom_line(aes(group=clone))+theme_bw()+scale_y_log10()


data <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1502/AF_thr_gained', sep= "\t", header=F, stringsAsFactors = F)
colnames(data) <- c('file', 'clone', 'MR')
remove <- c('platypus_nobin_00/all.MR_ov', 'platypus_nobin_cn2/all.MR_ov', 'platypus_nobin_indels/all.MR_ov', 'platypus_nobin/all.MR_ov')
data <- data[!data$file %in% remove,]

data$thr <- gsub('platypus_nobin_',  '', data$file)
data$thr <- gsub('/all.MR_ov',  '', data$thr)
data[data$thr=="cn3", 'thr'] <- 0
data[data$thr=="003", 'thr'] <- 0.03
data[data$thr=="006", 'thr'] <- 0.06
data[data$thr=="01", 'thr'] <- 0.1
data[data$thr=="03", 'thr'] <- 0.3
data[data$thr=="05", 'thr'] <- 0.5
data[data$thr=="07", 'thr'] <- 0.7
data[data$thr=="09", 'thr'] <- 0.9


ggplot(data=data, aes(y=MR, color=clone, x=thr))+geom_point()+geom_line(aes(group=clone))+theme_bw()+ylab('# gained muts')


#egrassi@godot:/scratch/trcanmed/AF_spectra/dataset/CRC1502$ for f in platypus_nobin_0*/*-1-*calls.tsv.gz; do zcat $f | sed 1d | bawk '$14==1' | wc -l | bawk -v n=$f '{print n,$0}'; done  > called_clone_1
#egrassi@godot:/scratch/trcanmed/AF_spectra/dataset/CRC1502$ for f in platypus_nobin_0*/*-0.calls.tsv.gz; do zcat $f | sed 1d | bawk '$14==1' | wc -l | bawk -v n=$f '{print n,$0}'; done  > called_clone_0

data <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1502/called_clone_1', sep= "\t", header=F, stringsAsFactors = F)
colnames(data) <- c('file', 'MR')
dd <- strsplit(data$file, '/')
data$file <- sapply(dd, function(x) {x[[1]]})
data$clone <- sapply(dd, function(x) {x[[2]]})


remove <- c('platypus_nobin_00', 'platypus_nobin_cn2', 'platypus_nobin_indels', 'platypus_nobin')
data <- data[!data$file %in% remove,]


data$clone <- gsub('.calls.tsv.gz',  '', data$clone, fixed=TRUE)
data$thr <- gsub('platypus_nobin_',  '', data$file)
data[data$thr=="cn3", 'thr'] <- 0
data[data$thr=="003", 'thr'] <- 0.03
data[data$thr=="006", 'thr'] <- 0.06
data[data$thr=="01", 'thr'] <- 0.1
data[data$thr=="03", 'thr'] <- 0.3
data[data$thr=="05", 'thr'] <- 0.5
data[data$thr=="07", 'thr'] <- 0.7
data[data$thr=="09", 'thr'] <- 0.9


ggplot(data=data, aes(y=MR, color=clone, x=thr))+geom_point()+geom_line(aes(group=clone))+theme_bw()+ylab('# muts over thr')


data <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1502/called_clone_0', sep= "\t", header=F, stringsAsFactors = F)
colnames(data) <- c('file', 'MR')
dd <- strsplit(data$file, '/')
data$file <- sapply(dd, function(x) {x[[1]]})
data$clone <- sapply(dd, function(x) {x[[2]]})


remove <- c('platypus_nobin_00', 'platypus_nobin_cn2', 'platypus_nobin_indels', 'platypus_nobin')
data <- data[!data$file %in% remove,]


data$clone <- gsub('.calls.tsv.gz',  '', data$clone, fixed=TRUE)
data$thr <- gsub('platypus_nobin_',  '', data$file)
data[data$thr=="cn3", 'thr'] <- 0
data[data$thr=="003", 'thr'] <- 0.03
data[data$thr=="006", 'thr'] <- 0.06
data[data$thr=="01", 'thr'] <- 0.1
data[data$thr=="03", 'thr'] <- 0.3
data[data$thr=="05", 'thr'] <- 0.5
data[data$thr=="07", 'thr'] <- 0.7
data[data$thr=="09", 'thr'] <- 0.9


ggplot(data=data, aes(y=MR, color=clone, x=thr))+geom_point()+geom_line(aes(group=clone))+theme_bw()+ylab('# muts over thr')