library(pheatmap)
library(RColorBrewer)
setwd('/scratch/trcanmed/AF_spectra/')
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_2a_cosmic.pdf.Rdata')
bd <- data
annot_rows_bulk <- annot_rows
keep <- apply(bd, 1, function(x) {any(x>0.1)})
keep2 <- apply(bd, 1, function(x) {any(x >= 0.1)})

load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_2abis_cosmic.pdf.Rdata')
ma <- data
annot_rows_MA <- annot_rows

keep <- keep | apply(ma, 1, function(x) {any(x>0.1)})
keep2 <- keep2 | apply(ma, 1, function(x) {any(x >= 0.1)})
all(keep2 == keep)
bd <- bd[keep,]
ma <- ma[keep,]

pdf('mfag1.pdf', family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
pheatmap(bd, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows_bulk, annotation_colors = annot_colors,  
         color=brewer.pal(5, 'Blues'), breaks=seq(0, 0.6, by=0.1),
         cluster_rows=FALSE, annotation_legend=FALSE)
dev.off()

pdf('mfag2.pdf', family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
pheatmap(ma, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows_MA, annotation_colors = annot_colors,  
         color=brewer.pal(5, 'Blues'), breaks=seq(0, 0.6, by=0.1),
         cluster_rows=FALSE, annotation_legend=FALSE)
dev.off()
#### progress

library(pheatmap)
library(RColorBrewer)
setwd('/scratch/trcanmed/AF_spectra/')
load('/scratch/trcanmed/AF_spectra/dataset_fIANG/fig_2a_cosmic.pdf.Rdata')
bd <- data
annot_rows_bulk <- annot_rows
keep <- apply(bd, 1, function(x) {any(x>0.1)})
keep2 <- apply(bd, 1, function(x) {any(x >= 0.1)})

load('/scratch/trcanmed/AF_spectra/dataset_fIANG/fig_2abis_cosmic.pdf.Rdata')
ma <- data
annot_rows_MA <- annot_rows

keep <- keep | apply(ma, 1, function(x) {any(x>0.1)})
keep2 <- keep2 | apply(ma, 1, function(x) {any(x >= 0.1)})
all(keep2 == keep)
bd <- bd[keep,]
ma <- ma[keep,]

pdf('/scratch/trcanmed/AF_spectra/temp/sign_bulk.pdf', family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
pheatmap(bd, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows_bulk, annotation_colors = annot_colors,  
         color=brewer.pal(5, 'Blues'), breaks=seq(0, 0.6, by=0.1),
         cluster_rows=FALSE, annotation_legend=FALSE)
dev.off()

pdf('/scratch/trcanmed/AF_spectra/temp/sign_MA.pdf', family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
pheatmap(ma, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows_MA, annotation_colors = annot_colors,  
         color=brewer.pal(5, 'Blues'), breaks=seq(0, 0.6, by=0.1),
         cluster_rows=FALSE, annotation_legend=FALSE)
dev.off()