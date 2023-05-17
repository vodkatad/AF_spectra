library(ggplot2)
#d <- read.table('/scratch/trcanmed/AF_spectra/dataset/dnds_150_manual.tsv', sep="\t", header=FALSE, stringsAsFactors = FALSE)
d <- read.table('/scratch/trcanmed/AF_spectra/dataset/dnds_15010_manual.tsv', sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(d) <- c('name', 'x', 'x1', 'mean', 'lower', 'upper')
d$mut <- ifelse(grepl('0.12', d$name), 'subclonal', 'all')
d$model <- sapply(strsplit(d$name, "_"), function(x) {x[[5]][1]})
d$model <- gsub('.subclonal', '',  d$model, fixed=TRUE)
d$model_mut <- paste0(d$model, d$mut)
ggplot(data=d, aes(x=model_mut, color=mut, y=mean))+geom_point()+geom_errorbar(aes(x=model_mut, ymin=lower, ymax=upper))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1)


### WIP edu
#egrassi@godot:/scratch/trcanmed/AF_spectra/local/share/data$ cut -f 1,2,3,4,5,6 MARCATUREFILE_ripulito.txt > MARCATUREFILE_ripulitoo.txt
library(stringr)
d <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/MARCATUREFILE_ripulitoo.txt', sep="\t", header=TRUE, stringsAsFactors = FALSE)
d$date <-  as.Date(d$DATA, format =  "%Y/%m/%d")
d$smodel <- str_match(d$CASO, "(CRC\\d+)[\\sA-Z]?")[,2]
d <- d[order(d$date),]
d2 <- d[d$smodel == "CRC0282",]
#ggplot(data=d, aes(x=date, y=EDU.FRA, group=smodel, color=smodel))+geom_point()+geom_line()+theme_bw()


ggplot(data=d2, aes(x=date, y=EDU.FRA))+geom_point()+geom_line()+theme_bw()


for (s in unique(d$smodel)) {
  if (!is.na(s)) {
    d2 <- d[d$smodel == s,]
    print(ggplot(data=d2, aes(x=date, y=EDU.FRA))+geom_point()+geom_line()+theme_bw()+ggtitle(s))
  }
}