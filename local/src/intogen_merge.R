
muts_f <- snakemake@input[["muts"]]
drivers_f <- snakemake@input[["drivers"]]
out_f <- snakemake@output[[1]]

muts <- read.table(muts_f, sep='\t', header=T, stringsAsFactors = F)
drivers <- read.table(drivers_f, sep='\t', header=T, stringsAsFactors = F)
#drivers <- drivers[, c('SYMBOL', 'CANCER_TYPE', 'ROLE')]
drivers <- drivers[, c('SYMBOL', 'CANCER_TYPE')]

## collapse
drivers <- drivers[drivers$SYMBOL %in% muts$Gene.refGene, ]

if (nrow(drivers) != 0) {
    collapsecol <- function(x, col, data) {
        data <- data[data$SYMBOL == x, ]
        return(c(x, paste0(data[,col], collapse=';')))
    }

    udrivers <- as.data.frame(t(sapply(unique(drivers$SYMBOL), collapsecol, 'CANCER_TYPE', drivers)))

    colnames(udrivers) <- c('SYMBOL', 'CANCER_TYPE')
} else {
    udrivers <- data.frame(SYMBOL=c('NOT_A_GENE'), CANCER_TYPE=c('')) #vabbe
}
save.image('meh.Rdata')
res <- merge(muts, udrivers, by.x='Gene.refGene', by.y='SYMBOL', all.x=T)
#res[!is.na(res$CANCER_TYPE), 'CANCER_TYPE'] <- ''
#res[!is.na(res$ROLE), 'ROLE'] <- ''

res$Intogen_driver <- ifelse(res$CANCER_TYPE != '', 'yes', 'no')
colnames(res)[colnames(res)=='CANCER_TYPE'] <- "Intogen_Cancer_Type"
#colnames(res)[colnames(res)=='ROLE'] <- "Intogen_Role"

write.table(res, file=out_f, sep='\t', quote=F, row.names=FALSE)

