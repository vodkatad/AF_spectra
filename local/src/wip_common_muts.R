# nonsyn_2nd_f <- '/scratch/trcanmed/AF_spectra/datasetV2/CRC1502_clones_all/platypus_nobin_00/nonsyn.binary.tsv.gz'
# nonsyn_1st_f <- '/scratch/trcanmed/AF_spectra/datasetV2/CRC1502/platypus_nobin_00/nonsyn.binary.tsv.gz'
# muts_1st_f <- '/scratch/trcanmed/AF_spectra/datasetV2/CRC1502/platypus_nobin_00/CRC1502-09.annovar.gz'
# muts_2nd_f <- '/scratch/trcanmed/AF_spectra/datasetV2/CRC1502_clones_all/platypus_nobin_00/CRC1502-09.annovar.gz'
# 
# target <- 'CRC1502-09C'
# targetT1 <- 'CRC1502-09-1-C'
# nbrother <- 3
# 
# outfile <- '/scratch/trcanmed/AF_spectra/datasetV2/CRC1502-09C_CRC1502-09-1-C.snv.private.tsv'
nonsyn_2nd_f  <- snakemake@input[['second_genes']]
nonsyn_1st_f  <- snakemake@input[['first_genes']]
muts_1st_f  <- snakemake@input[['first_muts']]
muts_2nd_f  <- snakemake@input[['second_muts']]

target <- snakemake@params[['target']]
targetT1 <- snakemake@params[['targetT1']]
nbrother <- snakemake@params[['nbrothers']]

outfile <- snakemake@output[['muts']]

second <- read.table(gzfile(nonsyn_2nd_f), sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)
first <- read.table(gzfile(nonsyn_1st_f), sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)


target_genes_mat <- second[grepl(target, rownames(second)),]
second_not_target <- second[!grepl(target, rownames(second)),]

shared_genes <- colnames(target_genes_mat)[colSums(target_genes_mat)==3]

mutated_in_any <- unique(c(colnames(first)[colSums(first)>0],colnames(second_not_target)[colSums(second_not_target)>0]))
genes_private_clones <- setdiff(shared_genes, mutated_in_any)

t2_private <- read.table(gzfile(muts_2nd_f), sep="\t", header=FALSE, stringsAsFactors = F)

nonsyn <- c('nonsynonymous SNV',
            'stopgain')
t2_private <- t2_private[t2_private$V8 %in% genes_private_clones & t2_private$V7 == 'exonic' & t2_private$V10 %in% nonsyn ,] #& t2_private$V1 == targetT1,]

if (length(unique(t2_private$V8))!=length(genes_private_clones)) {
  save.image('llama.Rdata')
  stop('assumptions not met regarding genes-muts private')
}
if (nrow(t2_private) != length(genes_private_clones)*nbrother) {
  save.image('llama.Rdata')
  stop('assumptions not met regarding genes-muts private')
}
# Now we need the ones found in the three _and_ their respective T1 but not in others.
# gained in T2 clearly do not include this!

# Gained in T1
t1 <- first[rownames(first)==targetT1, , drop=FALSE]
othert1 <- first[rownames(first)!=targetT1, , drop=FALSE]

genes_t1 <- colnames(t1)[t1[1,]!=0]
genes_othert1 <- colnames(othert1)[colSums(othert1)!=0]

genes_t1 <- setdiff(genes_t1, genes_othert1)
# go back to single muts
t1_all <- read.table(gzfile(muts_1st_f), sep="\t", header=FALSE, stringsAsFactors = F)

nonsyn <- c('nonsynonymous SNV',
            'stopgain')
t1_all <- t1_all[t1_all$V8 %in% genes_t1 & t1_all$V7 == 'exonic' & t1_all$V10 %in% nonsyn ,] #& t1_all$V1 == targetT1,]
# exonic is not in the rule
# limit on targetT1 not needed considering that we extracted the gained _only_ in our targetT1


if (nrow(t1_all) != length(genes_t1)) {
  save.image('llama.Rdata')
  stop('assumptions not met regarding genes-muts T1')
}

if (nrow(t1_all) > 0) {
  t1_all$status <- 'inherited'
}
if (nrow(t2_private) > 0) {
  t2_private$status <- 'denovo_private'
}
res <- rbind(t1_all, t2_private)
if (ncol(res) == 15) {
  colnames(res) <- c('clone', 'chr', 'b', 'e', 'region', 'ref', 'alt','gene', 'pass','effect','protein','dbSNP', 'cosmic', 'pass2','status')
}
res$clone <- NULL
res <- unique(res)
write.table(res, file=outfile, quote=FALSE, sep="\t", row.names = FALSE)
