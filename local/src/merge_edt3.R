library(reshape)

SNV_f  <- snakemake@input[['SNV']]
indel_f  <- snakemake@input[['indel']]
glen_f  <- snakemake@input[['glen']]
gens_f  <- snakemake@input[['gens']]
SNVn_f  <- snakemake@input[['SNVn']]
indeln_f  <- snakemake@input[['indeln']]


outtsv_f <- snakemake@output[['tsv']]

save.image(paste0(outtsv_f, '.Rdata'))

SNV <- read.table(SNV_f, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(SNV) <- c('id', 'MR_SNVs')

indel <- read.table(indel_f, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(indel) <- c('id', 'MR_indels')

SNVn <- read.table(SNVn_f, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(SNVn) <- c('id', 'N_SNVs')

indeln <- read.table(indeln_f, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(indeln) <- c('id', 'N_indels')


glen <- read.table(glen_f, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(glen) <- c('id', 'Genome_Length')

gens <- read.table(gens_f, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(gens) <- c('model', 'id', 'n_generations', 'measure', 'days')

## We need to get rid of 2n round QC failed CRC0282 from indels because their conf its conf still listed them.
## Keep everything there is in SNV, our main reference
n <- nrow(SNV)

m <- merge(SNV, indel, by="id")
mn <- merge(SNVn, indeln, by="id")
m <- merge(m, mn, by="id")
our <- merge(m, glen, by="id")

n2 <- nrow(our)
stopifnot(n2==n)

our$model <- sapply(our$id, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
our$clone <- sapply(our$id, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
our$clone2 <- sapply(our$id, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
our$model_clone <- paste0(our$model, "_", our$clone)

gens$model_clone <- paste0(gens$model, "_", gens$id)
wide_gens <- cast(gens, value="n_generations", formula=as.formula("model_clone~measure"))
colnames(wide_gens) <- c('model_clone', 'Generations_cell_doublings', 'Generations_EDU')
m3 <- merge(our, wide_gens, by="model_clone")

n3 <- nrow(m3)
stopifnot(n3==n)
# 25 (T2 clones) + 73 = 98 yay

m3$model_clone <- NULL
m3$clone <- NULL
m3$clone2 <- NULL

m3$model <- paste0(m3$model, ifelse(!grepl('\\d$', m3$model), '', ifelse(m3$model=="CRC0282", 'PR', 'LM')))

m3 <- m3[,c( 'model','id','MR_SNVs','MR_indels','N_SNVs','N_indels' ,'Genome_Length', 'Generations_EDU', 'Generations_cell_doublings')]

write.table(m3, file=outtsv_f, sep="\t", quote=FALSE, row.names = FALSE)
