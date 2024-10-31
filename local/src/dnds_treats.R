library(ggplot2)

d <- read.table('/scratch/trcanmed/AF_spectra/dataset/all_dnds_n.tsv', sep="\t", stringsAsFactors=FALSE, header=F)
colnames(d) <- c('class', 'n', 'name')

d$exp <- sapply(d$name, function(x) {strsplit(x, '/')[[1]][6]})
d$exp <- gsub('_platypus_nobin', '', d$exp)

tot <- sapply(unique(d$exp), function(x){dd <- d[d$exp==x,]; sum(dd$n)})
app <- data.frame(exp=unique(d$exp), tot=tot)
ggplot()+geom_col(data=d, aes(y=n, fill=class, x=exp), position='fill')+theme_bw(base_size=20)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) + geom_text(data=app, aes(x=exp, y=1.03,label=tot))
