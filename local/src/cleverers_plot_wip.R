files <- c('/scratch/trcanmed/AF_spectra/dataset/cleverers/CRC2566/all_cosmic_fit.tsv',
           '/scratch/trcanmed/AF_spectra/dataset/cleverers/CRC2573/all_cosmic_fit.tsv',
           '/scratch/trcanmed/AF_spectra/dataset/cleverers/CRC2608/all_cosmic_fit.tsv')


load_polish <- function(filename) {
  d <- read.table(filename, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)
  print(d[nrow(d),])
  d <- d[-nrow(d),]
  smodel <- strsplit(filename, '/')[[1]][7]
  colnames(d) <- paste0(colnames(d), "_", smodel)
  d
}

alld <- sapply(files, load_polish)
names(alld) <- NULL
cdata <- do.call(cbind, alld)

data <- t(as.matrix(cdata))
annot_rows <- data.frame(row.names=rownames(data))

annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'_', fixed=T), function(x){ x[2] }))
#annot_rows$clone <- as.factor(unlist(lapply(strsplit(rownames(annot_rows),'.', fixed=T), function(x){ paste0(x[2], '.', x[3]) })))

#col <- brewer.pal(3,'Dark2')
#names(col) <- levels(annot_rows$sample)
# cosmic specific XXX WARNING

colnames(data) <- unlist(lapply(strsplit(colnames(data),'.', fixed=TRUE), function(x){ x[length(x)] }))

#cbPalette2 <- unlist(strsplit(palette, ','))

#names(cbPalette2) <- levels(annot_rows$model)
#annot_colors <- list(sample=col, model=cbPalette2)
pheatmap(t(data), fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE,  annotation_col=annot_rows,   color=brewer.pal(9, 'PuBu'), width=6, height=8.3)#, filename="clevers.pdf")
