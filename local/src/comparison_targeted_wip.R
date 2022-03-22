library(ggplot2)

#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect$ bawk '$6!=0' merged_longformat_wtiers.tsv | grep CRC0282 | cut -f 2,6 | tr ":" "\t"  > ~/282_biobanca.tsv 
#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect$ bawk '$6!=0' merged_longformat_alltiers.tsv | grep CRC0282 | cut -f 2,6 | tr ":" "\t"  > ~/282_biobanca_all.tsv 
#egrassi@godot:/scratch/trcanmed/AF_spectra/local/share/data/CRC0282$ zcat CRC0282-01-1-A.pass.vcf.gz | grep -v "^#" | filter_1col 2 <(cut -f 2 ~/282_biobanca.tsv) | cut -f 1,2,4,5,10,11

#  egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect$ wc -l ~/282_alltiers_check.tsv
#35 /home/egrassi/282_alltiers_check.tsv
#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect$ wc -l ~/282_wtiers_check.tsv
# 18 ..

#targ <- read.table('~/282_biobanca.tsv', sep="\t", header=F, stringsAsFactors = F)
#MA <- read.table('~/282_wtiers_check.tsv', sep="\t", header=F, stringsAsFactors = F)


targ <- read.table('~/282_biobanca_all.tsv', sep="\t", header=F, stringsAsFactors = F)
MA <- read.table('~/282_alltiers_check.tsv', sep="\t", header=F, stringsAsFactors = F)

colnames(targ) <- c('chr', 'pos', 'ref', 'alt', 'af_targ')
rownames(targ) <- paste0(targ$chr, ":", targ$pos, ":", targ$ref, ":", targ$alt)
colnames(MA) <- c('chr', 'pos', 'ref', 'alt', 'clone', 'matchnorm')
rownames(MA) <- paste0(MA$chr, ":", MA$pos, ":", MA$ref, ":", MA$alt)
MA$af_MA <- as.numeric(sapply(strsplit(MA$clone, ":"), '[[', 3))

m <- merge(targ, MA, all=TRUE)
m[is.na(m$af_targ),'af_targ'] <- 0
m[is.na(m$af_MA),'af_MA'] <- 0

ggplot(data=m, aes(af_MA, af_targ))+geom_point()+theme_bw()+theme(text=element_text(size=20))
