acc_f  <- snakemake@input[['acc']]
essential_f  <- snakemake@input[['essential']]

log_f <- snakemake@log[['log']]
outtsv_f <- snakemake@output[['outtsv']]


g_ess <- read.table(essential_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(g_ess) <- 'gene_symbol'

res <- read.table(acc_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)

colnames(res) <- gsub('.', "-", colnames(res), fixed=TRUE)
res[is.na(res)] <- ''

find_essential <- function(vec, gs) {
  res <- rep(FALSE, length(vec))
  for (i in seq(1, length(vec))) {
    genes <- unlist(strsplit(vec[i], ';', fixed=TRUE))
    res[i] <- any(genes %in% gs)
  }
  res
}

res$mutid <- ifelse(find_essential(res$gene, g_ess$gene_symbol), 'essential', 'non-essential')

res <- res[order(res$chr, res$b),]
colnames(res)[seq(1, 9)] <- c('Chromosome', 'Coordinate hg38', 'Ref seq', 'Alt seq', 'Genomic context', 'Gene', 
                  'Effect', 'Protein change', 'Essential gene')

write.table(res, file=outtsv_f, quote=FALSE, sep="\t", row.names = FALSE)

sink(log_f)
print(nrow(res))
sink()
