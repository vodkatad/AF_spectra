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
annot_rows$mut <- ifelse(unlist(lapply(strsplit(rownames(annot_rows),'_', fixed=T), function(x){ x[1] }))!="common", 'Private', 'Shared')

#annot_rows$clone <- as.factor(unlist(lapply(strsplit(rownames(annot_rows),'.', fixed=T), function(x){ paste0(x[2], '.', x[3]) })))

#col <- brewer.pal(3,'Dark2')
#names(col) <- levels(annot_rows$sample)
# cosmic specific XXX WARNING

colnames(data) <- unlist(lapply(strsplit(colnames(data),'.', fixed=TRUE), function(x){ x[length(x)] }))

#cbPalette2 <- unlist(strsplit(palette, ','))

#names(cbPalette2) <- levels(annot_rows$model)
#annot_colors <- list(sample=col, model=cbPalette2)

pheatmap(t(data), fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE,  annotation_col=annot_rows,   color=brewer.pal(9, 'PuBu'), width=4.13, height=5.7)


pheatmap(t(data), fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE,  annotation_col=annot_rows,   color=brewer.pal(9, 'PuBu'), width=4.13, height=5.7, filename="cleverers_annoted.pdf")


pc <- 0.0001
ratiodf <- data.frame(row.names=rownames(data), ratio = log10((data[,'8']+pc)/(data[,'1']+pc)))

ratiodf$mut <- ifelse(unlist(lapply(strsplit(rownames(ratiodf),'_', fixed=T), function(x){ x[1] }))!="common", 'Private', 'Shared')
ratiodf$model <- unlist(lapply(strsplit(rownames(ratiodf),'_', fixed=T), function(x){ x[2] }))
ratiodf$id <- rownames(ratiodf)



ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

ggplot(data=ratiodf, aes(x=reorder(id, as.numeric(as.factor(model))), y=ratio, fill=mut))+geom_col(position="dodge")+
  theme_bw()+ctheme+
  xlab('Sample')+ylab("Log10 ratio") + ggtitle("Signature 8 / Signature 1")+
  scale_fill_manual(values=c('darkgoldenrod', 'darkgreen'))

ggsave('cleverers_ratio.pdf', width=8, height=5, units="in")
