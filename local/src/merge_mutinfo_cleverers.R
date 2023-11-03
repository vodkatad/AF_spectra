ovcnok_f  <- snakemake@input[['ovcnok']]
nr_f  <- snakemake@input[['allnr']]

outtsv_f <- snakemake@output[['outtsv']]

common_private <- read.table(ovcnok_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(common_private) <- c('mutid', 'sample')

nr <- read.table(nr_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(nr) <- c('chr', 'b', 'e', 'mutid', 'details')

spli <- strsplit(nr$details, ':', fixed=TRUE)
nr$alt <- as.numeric(sapply(spli, function(x) {x[[3]][1]}))
nr$tot <- as.numeric(sapply(spli, function(x) {x[[2]][1]}))

nr$chr <- NULL
nr$b <- NULL
nr$e <- NULL
nr$details <- NULL
#nr$alt <- NULL
#nr$tot <- NULL

m <- merge(common_private, nr, by="mutid")
#CRC1307-02-1-A_CRC1307-02-0     chr1:14898042:A:G       26/74
m$sample <- paste0(m$sample, "_plh")
res <- m[, c('sample', 'mutid', 'alt', 'tot')]

# we need to put together common info
cres <- res[res$sample == "common_plh",]
pres <- res[res$sample != "common_plh",]

# even if there are strange things with the same nr/nv shared by all
# the average that we compute here should be right
avg_nr <- function(mutid, data) {
	subset <- data[data$mutid==mutid,]
	return(c('common_plh', mutid, mean(subset$alt), mean(subset$tot)))
}
avg <- sapply(unique(cres$mutid), avg_nr, cres)
commonavg <- as.data.frame(t(avg))
colnames(commonavg) <- colnames(pres) 
res <- rbind(pres, commonavg)
res$info <- paste0(res$alt, '/', res$tot)
res$alt <- NULL
res$tot <- NULL
write.table(res, file=outtsv_f, quote=FALSE, sep="\t", row.names = FALSE, col.names=FALSE)


