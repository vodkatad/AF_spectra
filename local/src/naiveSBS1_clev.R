#!/usr/bin/env Rscript
library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg19'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
outputtsv <- args[2]
outputbarpl <- args[3]
palette <- args[4]
log_f <- args[5]
rdata_f <- args[6]
outputperc <- args[7]

# infile is tsv vfc / samplename_bulk|vivo|vitro
vcfs <- read.table(input, sep="\t", header=FALSE, stringsAsFactors = FALSE)
vcf_files <- vcfs$V1
sample_names <- vcfs$V2

# I wanted to use expand in the rule that generates signinput_..., therefore single samples
# will create empty vcf if they do not have any in vitro/in vivo clone and we need to get rid of them here
keep_noempty <- file.size(vcf_files) != 0L
vcf_files <- vcf_files[keep_noempty]
sample_names <- sample_names[keep_noempty]

sink(log_f)
print("N samples:")
print(length(vcf_files))
table((unlist(lapply(strsplit(sample_names,'_'), function(x){ x[1] }))))
print("2nd round")
print(table(grepl('2nd_', sample_names)))
sink()

vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
mut_mat <- mut_mat[, !grepl('2nd', colnames(mut_mat))]
save.image(rdata_f)

wanted <- c('A[C>T]G','C[C>T]G','G[C>T]G','T[C>T]G')

totMuts <- as.data.frame(colSums(mut_mat))
mut_mat_SBS1 <- mut_mat[rownames(mut_mat) %in% wanted,]

totSbs1 <- as.data.frame(colSums(mut_mat_SBS1))

colnames(totSbs1) <- 'nSBS1'
colnames(totMuts) <- 'tot'
pdata <- merge(totSbs1, totMuts, by="row.names")
pdata$frac <- pdata$nSBS1 / pdata$tot

rownames(pdata) <- pdata$Row.names
pdata$Row.names <- NULL

pdata$sample <-  as.factor(unlist(lapply(strsplit(rownames(pdata),'_'), function(x){ x[length(x)] })))
pdata$model <- unlist(lapply(strsplit(rownames(pdata),'-'), function(x){ x[1] }))
pdata$model <- as.factor(unlist(lapply(strsplit(pdata$model,'_'), function(x){ x[1] })))

col <- brewer.pal(3,'Dark2')
names(col) <- levels(pdata$sample)

cbPalette2 <- unlist(strsplit(palette, ','))

names(cbPalette2) <- levels(pdata$model)
annot_colors <- list(sample=col, model=cbPalette2)

ggplot(data=pdata, aes(x=model,y=nSBS1, fill=model))+geom_col()+theme_bw()+facet_wrap(~sample)+
  geom_text(aes(label = nSBS1), hjust=0.5, angle=90, size=5)+
  scale_fill_manual(values=annot_colors$model)+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

ggsave(outputbarpl, width = 34, height=17, units="cm")

ggplot(data=pdata, aes(x=model,y=frac, fill=model))+geom_col()+theme_bw()+facet_wrap(~sample)+
  scale_fill_manual(values=annot_colors$model)+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

ggsave(outputperc, width = 34, height=17, units="cm")

write.table(pdata, file=outputtsv, sep="\t", quote=FALSE)

