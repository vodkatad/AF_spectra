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


cosmic <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "") # ???

setwd('/scratch/trcanmed/AF_spectra/dataset/')

#files <- c("CRC0327-02.signin.vcf","CRC0327-04.signin.vcf","CRC0327-08.signin.vcf","CRC0441-01.signin.vcf","CRC0441-03.signin.vcf","CRC0441-10.signin.vcf","CRC1078-02.signin.vcf","CRC1078-07.signin.vcf","CRC1078-09.signin.vcf","CRC1307-02.signin.vcf","CRC1307-08.signin.vcf","CRC1307-09.signin.vcf","CRC1502-03.signin.vcf","CRC1502-08.signin.vcf","CRC1502-09.signin.vcf","CRC1502-10.signin.vcf","CRC1599LM-01.signin.vcf","CRC1599LM-03.signin.vcf","CRC1599LM-07.signin.vcf","CRC1599PR-01.signin.vcf","CRC1599PR-10.signin.vcf")
#sample_names <- c("CRC0327-02","CRC0327-04","CRC0327-08","CRC0441-01","CRC0441-03","CRC0441-10","CRC1078-02","CRC1078-07","CRC1078-09","CRC1307-02","CRC1307-08","CRC1307-09","CRC1502-03","CRC1502-08","CRC1502-09","CRC1502-10","CRC1599LM-01","CRC1599LM-03","CRC1599LM-07","CRC1599PR-01","CRC1599PR-10")
input <- 'mehh'

vcfs <- read.table(input, sep="\t", header=FALSE, stringsAsFactors = FALSE)
vcf_files <- vcfs$V1
sample_names <- vcfs$V2

vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)


nrun_estimate <- 50
nrun <- 150
seed <- 123456
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:5, method="brunet", nrun=nrun_estimate, seed=seed)
plot(estimate)

nsign <- 3
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)

base <- 'toplevel'
names_sign <- paste0("sign_", base, "_", seq(1, nsign))
colnames(nmf_res$signatures) <- names_sign
rownames(nmf_res$contribution) <- names_sign

plot_96_profile(nmf_res$signatures, condensed = TRUE)

mypc <- function (contribution, signatures, index = c(), coord_flip = FALSE, mode = "relative", palette = c()) 
{
  if (!(mode == "relative" | mode == "absolute")) 
    stop("mode parameter should be either 'relative' or 'absolute'")
  if (length(index > 0)) {
    contribution = contribution[, index]
  }
  Sample = NULL
  Contribution = NULL
  Signature = NULL
  if (mode == "relative") {
    m_contribution = melt(contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    plot = ggplot(m_contribution, aes(x = factor(Sample), 
                                      y = Contribution, fill = factor(Signature), order = Sample)) + 
      geom_bar(position = "fill", stat = "identity", colour = "black") + 
      labs(x = "", y = "Relative contribution") + theme_bw() + 
      theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) + 
      theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())+
      theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1))
  }
  else {
    if (missing(signatures)) 
      stop(paste("For contribution plotting in mode 'absolute':", 
                 "also provide signatures matrix"))
    total_signatures = colSums(signatures)
    abs_contribution = contribution * total_signatures
    m_contribution = melt(abs_contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    plot = ggplot(m_contribution, aes(x = factor(Sample), 
                                      y = Contribution, fill = factor(Signature), order = Sample)) + 
      geom_bar(stat = "identity", colour = "black") + labs(x = "", 
                                                           y = "Absolute contribution \n (no. mutations)") + 
      theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                         panel.grid.major.x = element_blank()) + theme(panel.grid.minor.y = element_blank(), 
                                                                       panel.grid.major.y = element_blank())+
      theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1))
  }
  if (length(palette) > 0) 
    plot = plot + scale_fill_manual(name = "Signature", values = palette)
  else plot = plot + scale_fill_discrete(name = "Signature")
  if (coord_flip) 
    plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
  else plot = plot + xlim(levels(factor(m_contribution$Sample)))
  return(plot)
}

mypc(nmf_res$contribution, nmf_res$signature, mode = "relative")

plot_contribution_heatmap(nmf_res$contribution,sig_order = names_sign,cluster_samples=F)

sp_url <- paste(cosmic, sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe> 
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names>
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average") # store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
plot_cosine_heatmap(cos_sim_samples_signatures,col_order = cosmic_order,cluster_rows = TRUE)

us <- nmf_res$signatures
colnames(us) <- names_sign
hclust_cosmic_us = cluster_signatures(cbind(us, cancer_signatures), method = "average") # store signatures in new order
png('hcust.png')
plot(hclust_cosmic_us)
graphics.off()
# put together vitro-vivo of each starting clone and redo (to get more power in reconstructing signatures?)


# 5 - 18 with k=2 and 1 sample for each starting clone - with k = 3 signature 8 pops out
# 5 - 18 - 8 with k=3 and 1 sample for each model (vivo/vitro separated)
