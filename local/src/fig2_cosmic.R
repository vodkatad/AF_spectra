sign_f  <- snakemake@input[['sign']]
order_f  <- snakemake@input[['order']]
colors <- snakemake@input[['palette']]
wanted <- snakemake@params[['wanted']]
legend <- snakemake@params[['legend']]
outplot <- snakemake@output[['plot']]
addlm <- snakemake@params[['addlm']]

legendbool <- FALSE
if (legend=="yes") {
  legendbool <- TRUE
}
#theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(RColorBrewer)
library(pheatmap)

#load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

data <- read.table(sign_f, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)

if (!grepl('2nd', wanted)) {
  data <- data[grepl(wanted, rownames(data)),]
  data <- data[!grepl('2nd', rownames(data)),]
} else {
  stop('Still to be implemented')
}
# now we have each model only once so we can user match to get the right order
# not for clevers :(
data$model <- unlist(lapply(strsplit(rownames(data),'-'), function(x){ x[1] }))
data$model <- as.character(unlist(lapply(strsplit(data$model,'_'), function(x){ x[1] })))
data$r <- rownames(data)
orderdf <- read.table(order_f, sep="\t", quote="", header=TRUE, stringsAsFactor=FALSE)
orderdf$n <- seq(1, nrow(orderdf))
m <- merge(data, orderdf, by="model")
m <- m[order(m$n),]
data <- m[, colnames(data)]
data$model <- NULL
rownames(data) <- m$r
data$r <- NULL

annot_rows <- data.frame(row.names=rownames(data))

#annot_rows$sample <-  as.factor(unlist(lapply(strsplit(rownames(annot_rows),'_'), function(x){ x[length(x)] })))
annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] }))
annot_rows$model <- as.factor(unlist(lapply(strsplit(annot_rows$model,'_'), function(x){ x[1] })))
# add back LMX P
if (addlm!= "no") {
  annot_rows$model <- paste0(annot_rows$model, ifelse(!grepl('\\d$', annot_rows$model), '', ifelse(annot_rows$model=="CRC0282", 'PR', 'LM')))
  names(pal) <- paste0(names(pal), ifelse(!grepl('\\d$', names(pal)), '', ifelse(names(pal)=="CRC0282", 'PR', 'LM')))
}
colnames(data) <- gsub('X', 'SBS', colnames(data))
data <- t(data)  

annot_colors <- list(model=pal)
pdf(outplot, family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
ph <- pheatmap(data, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  
         color=brewer.pal(9, 'Blues'), breaks=seq(0, 0.6, by=0.07),
         cluster_rows=FALSE, annotation_legend=FALSE)
ph$gtable[[1]][[1]]$children[[1]]$gp$lwd <- 0.001
ph
graphics.off()
if (legendbool) {
  pdf(paste0('legend_', outplot), family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
  pheatmap(data, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
                 cluster_cols=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  
                 color=brewer.pal(9, 'Blues'), breaks=seq(0, 0.6, by=0.07),
                 cluster_rows=FALSE, annotation_legend=TRUE)
  graphics.off()
  
}
#pdf('test.pdf', family='sans')
#pheatmap(data, fontsize_row = 1.5, fontsize_col=1.5, fontsize.number=1.5, show_colnames = TRUE,  cluster_cols=FALSE, 
# annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), cluster_rows=FALSE, legend=TRUE)
#graphics.off()

save.image(paste0(outplot, '.Rdata'))
