
library(pheatmap)

d <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset/SNV_nonsyn.binary.tsv.gz'), sep="\t",header=TRUE, row.names = 1)
genes_n <- colSums(d)

cancer <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/intogen_symbols', header=F, sep="\t")
essential <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/common_iorio_crisp_essential',header=F, sep="\t")

hist(genes_n, n = 15)
#pheatmap(t(d), show_rownames = FALSE, legend = FALSE)

annot <- data.frame(row.names=names(genes_n))
annot$cancer <- ifelse(rownames(annot) %in% cancer$V1, 'yes', 'no')
colors <- list(cancer = c(yes = "black", no="darkgoldenrod"))

pheatmap(t(d), show_rownames = FALSE, legend = FALSE, annotation_row=annot, annotation_colors=colors)


annot <- data.frame(row.names=names(genes_n))
annot$essential <- ifelse(rownames(annot) %in% essential$V1, 'yes', 'no')
colors <- list(essential = c(yes = "black", no="darkgoldenrod"))

pheatmap(t(d), show_rownames = FALSE, legend = FALSE, annotation_row=annot, annotation_colors=colors)

# trying order based on essentiality
td <- t(d)
annot <- annot[order(annot$essential),, drop=F]
td <- td[match(rownames(annot), rownames(td)),]
pheatmap(td, show_rownames = FALSE, cluster_rows=F,legend = FALSE, annotation_row=annot, annotation_colors=colors)

annot <- data.frame(row.names=names(genes_n))
annot$cancer <- ifelse(rownames(annot) %in% cancer$V1, 'yes', 'no')
annot <- annot[order(annot$cancer),, drop=F]
td <- td[match(rownames(annot), rownames(td)),]
colors <- list(cancer = c(yes = "black", no="darkgoldenrod"))
pheatmap(td, show_rownames = FALSE, cluster_rows=F,legend = FALSE, annotation_row=annot, annotation_colors=colors)

## CRC1502 investigations on clone 9
CRC1502 <- d[grepl('CRC1502',rownames(d), fixed=F), ]
rs <- colSums(CRC1502)
remove <- names(rs[rs==0])
CRC1502 <- CRC1502[, !colnames(CRC1502) %in% remove]
n <- rownames(CRC1502); n[grepl('CRC1502-09', n)]
clone9 <- CRC1502[grepl('CRC1502-09', n),]
others <- CRC1502[!grepl('CRC1502-09', n),]

list_mut <- function(data) {
  rs <- colSums(data)
  return(names(rs[rs>0]))
# }

list_all_mut <- function(data) {
  rs <- colSums(data)
  return(names(rs[rs==nrow(data)]))
}

CRC1502_base <- d[grepl('CRC1502-09-1-C',rownames(d), fixed=F), ]

all <- colSums(CRC1502)
aall <- colSums(d)
all[names(all) %in% list_all_mut(CRC1502_base)]
aall[names(aall) %in% list_all_mut(CRC1502_base)]

