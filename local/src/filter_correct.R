genes_f  <- snakemake@input[['genes']]
wilcox_f  <- snakemake@input[['wilcox']]

outfile <- snakemake@output[['tsv']]

genes <- read.table(genes_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
wilcox <- read.table(wilcox_f, sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)

wilcox$padj_strict <- NULL
wilcox$padj <- NULL

wilcox <- wilcox[rownames(wilcox) %in% genes$gene,]
wilcox$padj <- p.adjust(wilcox$pvalue, method="BH")

res <- cbind(data.frame(gene=rownames(wilcox)), wilcox)

write.table(res, file=outfile, sep="\t", row.names=FALSE, quote=FALSE)
