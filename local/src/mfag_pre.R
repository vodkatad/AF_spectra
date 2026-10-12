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
#### progress + GPS erc ###############################################################################################

library(pheatmap)
library(RColorBrewer)
setwd('/scratch/trcanmed/AF_spectra/')
load('/scratch/trcanmed/AF_spectra/dataset_fIANG/fig_2a_cosmic.pdf.Rdata')
bd <- data
annot_rows_bulk <- annot_rows
#keep <- apply(bd, 1, function(x) {any(x>0.1)})
#keep2 <- apply(bd, 1, function(x) {any(x >= 0.1)})

load('/scratch/trcanmed/AF_spectra/dataset_fIANG/fig_2abis_cosmic.pdf.Rdata')
ma <- data
annot_rows_MA <- annot_rows

annot_colors_orig <- annot_colors
#keep <- keep | apply(ma, 1, function(x) {any(x>0.1)})
#keep2 <- keep2 | apply(ma, 1, function(x) {any(x >= 0.1)})
#all(keep2 == keep)
#bd <- bd[keep,]
#ma <- ma[keep,]

#### fish back MSI ####
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_2a_cosmic.pdf.Rdata')
bd2 <- data
annot_rows_bulk2 <- annot_rows
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_2abis_cosmic.pdf.Rdata')
ma2 <- data
annot_rows_MA2 <- annot_rows
annot_colors2 <- annot_colors
####
ma <- cbind(ma, ma2[, 'CRC0282_vitroMA', drop=F])
bd <- cbind(bd, bd2[, 'CRC0282_bulk', drop=F])

colnames(ma) <- gsub('CRC0282_vitroMA', 'CRC0282PR_vitroMA', colnames(ma), fixed=T)
colnames(bd) <- gsub('CRC0282_bulk', 'CRC0282PR_bulk', colnames(bd), fixed=T)

keep <- apply(bd, 1, function(x) {any(x>0.1)})
keep2 <- apply(bd, 1, function(x) {any(x >= 0.1)})
keep <- keep | apply(ma, 1, function(x) {any(x>0.1)})
keep2 <- keep2 | apply(ma, 1, function(x) {any(x >= 0.1)})
all(keep2 == keep)
bd <- bd[keep,]
ma <- ma[keep,]

# order according to MR
mod <- c("CRC0282PR", "CRC1599PR","CRC2826PRO","CRC3023PRO","CRC0327","CRC0441","CRCUECHPRO","CRC1502","CRC1599LM","CRC1078","CRC1307")
modma <- paste0(mod, '_vitroMA')
modb <- paste0(mod, '_bulk')

dim(bd)
dim(ma)

ma <- ma[, modma]
bd <- bd[, modb]

dim(bd)
dim(ma)

annot_rows_MA <- rbind(annot_rows_MA, annot_rows_MA2['CRC0282_vitroMA',, drop=F])
annot_rows_bulk <- rbind(annot_rows_bulk, annot_rows_bulk2['CRC0282_bulk',, drop=F])
rownames(annot_rows_MA) <- gsub('CRC0282_vitroMA', 'CRC0282PR_vitroMA', rownames(annot_rows_MA), fixed=T)
rownames(annot_rows_bulk) <- gsub('CRC0282_bulk', 'CRC0282PR_bulk', rownames(annot_rows_bulk), fixed=T)

annot_colors_m <- annot_colors_orig
annot_colors_m$model <- c(annot_colors_m$model, head(annot_colors2$model, n=1))

pdf('/scratch/trcanmed/AF_spectra/temp/sign_bulk.pdf', family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
pheatmap(bd, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows_bulk, annotation_colors = annot_colors_m,  
         color=brewer.pal(5, 'Blues'), breaks=seq(0, 0.6, by=0.1),
         cluster_rows=FALSE, annotation_legend=FALSE)
dev.off()

pdf('/scratch/trcanmed/AF_spectra/temp/sign_MA.pdf', family="sans")#, width=2.2, height=1.4) # resize by hand cause otherwise it will be a mess
pheatmap(ma, cellwidth=5.67, cellheight=5.67, fontsize_row = 5, fontsize_col=5, fontsize.number=5, show_colnames = FALSE, show_rownames = TRUE,  
         cluster_cols=FALSE, annotation_col=annot_rows_MA, annotation_colors = annot_colors_m,  
         color=brewer.pal(5, 'Blues'), breaks=seq(0, 0.6, by=0.1),
         cluster_rows=FALSE, annotation_legend=FALSE)
dev.off()
























###### ?


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

# order according to MR
mod <- c("CRC0282", "CRC1599PR","CRC2826PRO","CRC3023PRO","CRC0327","CRC0441","CRCUECHPRO","CRC1502","CRC1599LM","CRC1078","CRC1307")
modma <- paste0(mod, '_vitroMA')
modb <- paste0(mod, '_bulk')

dim(bd)
dim(ma)

ma <- ma[, modma]
bd <- bd[, modb]

dim(bd)
dim(ma)

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