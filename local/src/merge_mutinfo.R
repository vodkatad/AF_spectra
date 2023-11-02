bulk_f  <- snakemake@input[['bulk']]
acc_f  <- snakemake@input[['acc']]
essential_f  <- snakemake@input[['essential']]

log_f <- snakemake@log[['log']]
outtsv_f <- snakemake@output[['outtsv']]

save.image(paste0(outtsv_f, '.Rdata'))

#load('/scratch/trcanmed/AF_spectra/datasetV2/edt_dir/CRC1307_edt2.tsv.Rdata')
#setwd('/scratch/trcanmed/AF_spectra/datasetV2')

g_ess <- read.table(essential_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(g_ess) <- 'gene_symbol'

bulk <- read.table(bulk_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
acc <- read.table(acc_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)

ff <- intersect(acc$mutid, bulk$mutid)

# only bulk
nclones <- ncol(acc)-9
bulkonly <- setdiff(bulk$mutid, acc$mutid)
res <- bulk[bulk$mutid %in% bulkonly, ]
add_acc_nm <- data.frame(matrix(rep('', length(bulkonly)*nclones), nrow=length(bulkonly)))
colnames(add_acc_nm) <- colnames(acc)[seq(10, ncol(acc))]
res <- cbind(res, add_acc_nm)

#only acc
acconly <- setdiff(acc$mutid, bulk$mutid)
only_acc <- acc[acc$mutid %in% acconly, ]
only_acc[, colnames(bulk)[10]] <- rep('', nrow(only_acc))
only_acc <- only_acc[, c(seq(1, 9), ncol(only_acc), seq(10, ncol(only_acc)-1))]

res <- rbind(res, only_acc)

# both, we get annot from bulk
bulk_both <- bulk[bulk$mutid %in% ff, ]
acc_both <- acc[acc$mutid %in% ff, seq(10, ncol(acc))]
res <- rbind(res, cbind(bulk_both, acc_both))

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
print(length(ff))
print(length(bulkonly))
print(length(acconly))
sink()
