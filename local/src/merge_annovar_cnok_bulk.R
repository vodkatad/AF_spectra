annovar_f  <- snakemake@input[['annovar']]
mutinfo_f  <- snakemake@input[['mutinfo']]

outtsv_f <- snakemake@output[['outtsv']]
kind <- snakemake@params[['kind']]
name <- snakemake@params[['name']]

save.image(paste0(outtsv_f, '.Rdata'))

library(reshape)

#load('/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/platypus_nobin_00/mutinfo.tsv.gz.Rdata')
#setwd('/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/platypus_nobin_00')
d_ann <- read.table(gzfile(annovar_f), sep="\t", header=TRUE, stringsAsFactors = FALSE)
d_ann <- d_ann[,c(1,2,3,4,5,6,7,8,10,11)]
colnames(d_ann) <- c('id', 'chr', 'b', 'e', 'ref', 'alt', 'location', 'gene', 'effect', 'aa_change')
d_ann$id <- rep(name, nrow(d_ann))
d_ann$aa_change <- sapply(strsplit(d_ann$aa_change, ":", fixed=TRUE), function(x) {l <- length(x); return(x[[l]][1]) })
d_ann$mutid <- paste0('chr', d_ann$chr, ":", d_ann$b, ":", d_ann$ref, ":", d_ann$alt)
all_muts <- read.table(mutinfo_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(all_muts) <- c('sample', 'mutid', 'VAF')

## check that we have all the clones we want
all_muts$id <- sapply(strsplit(all_muts$sample, "_", fixed=TRUE), function(x) {l <- length(x); return(x[[1]][1]) })

l1 <- length(unique(all_muts$id))
l2 <- length(unique(d_ann$id))
l3 <- length(intersect(unique(all_muts$id), unique(d_ann$id)))
stopifnot(l1==l2)
stopifnot(l2==l3)

# to wide format for d_ann
d_ann <- unique(d_ann)
d_ann$value <- 1
wide_anno <- cast(d_ann, formula = "chr + b + e + ref + alt + location + gene + effect + aa_change + mutid ~ id", value="value", fun.aggregate=sum)
#summary(unlist(wide_anno[, 11:ncol(wide_anno)])))

wide_muts <- cast(all_muts, formula="mutid~id", value="VAF")#, fill=)

n1 <- nrow(wide_muts)
#if (kind == "SNV") {
#  stopifnot(length(intersect(wide_muts$mutid, wide_anno$mutid))==n1)
#} else {
#  missedov <- setdiff(wide_muts$mutid, wide_anno$mutid)
#  write.table(data.frame(mut=missedov), file=paste0(outtsv_f, ".miss"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
#in bulk we have all vars with all CN in the input
  common <- intersect(wide_muts$mutid, wide_anno$mutid)
  wide_anno <- wide_anno[wide_anno$mutid %in% common,]
  wide_muts <- wide_muts[wide_muts$mutid %in% common,]
  # otherwise the later stopifnot throws an error if we have a data.frame/cast.df
  wide_muts <- as.data.frame(wide_muts)
  wide_anno <- as.data.frame(wide_anno)
#}
# order wide_muts like wide_anno
stopifnot(nrow(wide_muts)==nrow(wide_anno))
wide_muts <- wide_muts[match(wide_anno$mutid, wide_muts$mutid),]
row.names(wide_muts) <- wide_muts$mutid

wide_muts$mutid <- NULL
# check that mutational info corresponds
stopifnot(!any(is.na(wide_muts[wide_anno[, 11:ncol(wide_anno), drop=FALSE] != 0])))

# we can plug in ref/alt instead of 1/0

res <- cbind(wide_anno[,seq(1, 10)], wide_muts)

res$e <- NULL
write.table(res, file=gzfile(outtsv_f), sep="\t", quote=FALSE, row.names = FALSE)
