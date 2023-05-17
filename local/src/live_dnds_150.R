library(ggplot2)
d <- read.table('/scratch/trcanmed/AF_spectra/dataset/dnds_150_manual.tsv', sep="\t", header=FALSE, stringsAsFactors = FALSE)

colnames(d) <- c('name', 'x', 'x1', 'mean', 'lower', 'upper')
d$mut <- ifelse(grepl('0.12', d$name), 'subclonal', 'all')
d$model <- sapply(strsplit(d$name, "_"), function(x) {x[[1]][1]})
