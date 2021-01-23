#!/usr/bin/env Rscript
library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)

#cosmic <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "") # ???

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
outputdir <- args[2]
#nsign <- as.numeric(args[3])
setwd(outputdir)

### load data and define support var with sanples, etc
files <- unlist(strsplit(input, ','))
sample_names <- unlist(lapply(files, function(x) {y <- basename(x); strsplit(y,"_")[[1]][1]} ))
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
base_clone <- unique(sapply(sample_names, function(x) {y<-strsplit(x, '-')[[1]]; return(paste0(y[1],'-',y[2]))}))
base <- unique(sapply(sample_names, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])}))

### basic mut spectrum
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
pdf("spectrum.pdf")
plot_spectrum(type_occurrences, CT = TRUE)
graphics.off()
pdf("samples_spectrum.pdf")
plot_spectrum(type_occurrences, CT = TRUE, by=sample_names)
graphics.off()

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)


#garb <- sapply(base_clone, function(x) { 
  #pdf(paste0(x, ".96_profile.pdf"))
 # ggsave(file=paste0(x, ".96_profile.pdf"), plot=plot_96_profile(mut_mat[,grepl(x, colnames(mut_mat))]))
  #graphics.off()
#} )


### Transcriptional effects - strand bias

tx <- genes(eval(parse(text = ref_transcriptome)))
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, tx)

strand_counts <- strand_occurrences(mut_mat_s, by=sample_names)
strand_bias <- strand_bias_test(strand_counts)


pdf(paste0("strand_biases.pdf"))
gb <- lapply(base_clone,function(x) {plot_strand_bias(strand_bias[grepl(x,strand_bias$group),])})
do.call(grid.arrange, gb)
graphics.off()



save.image("mut_pat2.Rdata")

chromosomes <- seqnames(get(ref_genome))[1:22]
   gb <- lapply(sample_names,function(x) {
   ggsave(file=paste0(x,"_waterfall.png"), plot=plot_rainfall(vcf = vcfs[[x]], title=x, chromosomes=chromosomes))
})

