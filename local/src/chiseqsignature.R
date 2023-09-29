library(reshape)
library(ggplot2)
library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(RColorBrewer)
library(pheatmap)
# data is mut_mat from mut_pat_signatures_fit.R
load('/scratch/trcanmed/AF_spectra/datasetV2/vitrovivobulk_SBS1.Rdata')
### cosmic
cosmic <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "") # ???
sp_url <- paste(cosmic, sep = "")
ref_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), ref_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe> 
ref_signatures = ref_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names>
row.names(ref_signatures) = ref_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
ref_signatures = as.matrix(ref_signatures[,4:33])
##
ff <- fit_to_signatures(mut_mat, ref_signatures)

wanted <- c('1', '8', '18')
#### route 1 ############
data <- as.data.frame(t(ff$contribution))
data <- data/rowSums(data)
wdata <- data[,colnames(data) %in% paste0('Signature.', wanted)]

normalized <- wdata/rowSums(wdata)
####

##### route 0 ############
#data <- as.data.frame(ff$contribution)
#wdata <- data[rownames(data) %in% paste0('Signature.', wanted),]
#twdata <- t(wdata)
#normalized <- twdata/rowSums(twdata)
##################

bulk <- normalized[grepl('_bulk', rownames(normalized)),]
vitro <- normalized[grepl('_vitroMA', rownames(normalized)),]

vitrom <- as.data.frame(colMeans(vitro))
colnames(vitrom) <- 'vitroMA'
bulkm <- as.data.frame(colMeans(bulk))
colnames(bulkm) <- 'bulk'

pdata <- as.data.frame(t(cbind(vitrom, bulkm)))
pdata$sample <- rownames(pdata) 
mpdata <- melt(pdata)

ggplot(data=mpdata, aes(x=sample, y=value, fill=variable))+geom_col()+theme_bw()

##
data <- as.data.frame(t(ff$contribution))
data <- data/rowSums(data)
wdata <- data[,colnames(data) %in% paste0('Signature.', wanted)]
wdata$class <- sapply(strsplit(rownames(wdata), "_"), '[[', 2)
wdata$id <- rownames(wdata)
mp <- melt(wdata)
ggplot(data=mp, aes(y=value,x=id, fill=class))+geom_col()+facet_wrap(~variable)
ggplot(data=mp, aes(y=value,fill=variable, x=class))+geom_violin()+geom_jitter()+facet_wrap(~variable)

bu <- table(mp[mp$value > 0.1 & mp$class=="bulk", 'variable'])
vu <- table(mp[mp$value > 0.1 & mp$class=="vitroMA", 'variable'])

fid <- data.frame(bulk=bu, vitro=vu)
fid$bulk.Var1 <- NULL
rownames(fid) <- fid$vitro.Var1
fid$vitro.Var1 <- NULL
