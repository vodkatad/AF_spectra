sign_f  <- snakemake@input[['sign']]
colors <- snakemake@input[['palette']]
what <- snakemake@wildcards[['what']]
outplot <- snakemake@output[['plot']]

legendbool <- TRUE
#theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(RColorBrewer)
library(pheatmap)

#load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

data <- read.table(sign_f, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)

save.image(paste0(outplot, '.Rdata'))

if (!grepl('baseline', what)) {
  data <- data[!grepl('baseline', rownames(data)),]
}

#data$model <- unlist(lapply(strsplit(rownames(data),'-'), function(x){ x[1] }))
data$model <- as.character(unlist(lapply(strsplit(rownames(data),'_'), function(x){ x[1] })))
data$treat <- as.character(unlist(lapply(strsplit(rownames(data),'_'), function(x){ x[2] })))
data$r <- rownames(data)
data$faketreat <- gsub('Afatinib', 'zzz', data$treat)
data <- data[order(data$model, data$faketreat),]

annot_rows <- data.frame(row.names=rownames(data))

#annot_rows$sample <-  as.factor(unlist(lapply(strsplit(rownames(annot_rows),'_'), function(x){ x[length(x)] })))
annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] }))
annot_rows$model <- as.factor(unlist(lapply(strsplit(annot_rows$model,'_'), function(x){ x[1] })))
if (!grepl('baseline', what)) {
  annot_rows$treat <- factor(data$treat,  levels=c('NT', 'Afatinib'))
} else {
  annot_rows$treat <- factor(data$treat,  levels=c('baseline', 'NT', 'Afatinib'))
  
}

data$faketreat <- NULL
data$model <- NULL
data$treat <- NULL
data$r <- NULL
colnames(data) <- gsub('X', 'SBS', colnames(data))
data <- t(data)  

if (!grepl('baseline', what)) {
  treatcol <- c('#C0C0C0', '#FF8000')
  names(treatcol) <- c('NT',  'Afatinib')
} else{
  treatcol <- c('#C0C0C0', '#404040', '#FF8000')
  names(treatcol) <- c('NT', 'baseline', 'Afatinib')
}
annot_colors <- list(model=pal, treat=treatcol)
pdf(outplot, family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
ph <- pheatmap(data, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  
         color=brewer.pal(5, 'Blues'), breaks=seq(0, 0.4, by=0.1),
         cluster_rows=FALSE, annotation_legend=FALSE)
ph$gtable[[1]][[1]]$children[[1]]$gp$lwd <- 0.001
ph
graphics.off()
if (legendbool) {
  pdf(paste0('legend_', outplot), family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
  pheatmap(data, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
                 cluster_cols=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  
                 color=brewer.pal(5, 'Blues'), breaks=seq(0, 0.4, by=0.1),
                 cluster_rows=FALSE, annotation_legend=TRUE)
  graphics.off()
  
}

#> seq(0, 0.6, by=0.1)
#[1] 0.0 0.1 0.2 0.3 0.4 0.5 0.6

#pdf('test.pdf', family='sans')
#pheatmap(data, fontsize_row = 1.5, fontsize_col=1.5, fontsize.number=1.5, show_colnames = TRUE,  cluster_cols=FALSE, 
# annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), cluster_rows=FALSE, legend=TRUE)
#graphics.off()

save.image(paste0(outplot, '.Rdata'))
