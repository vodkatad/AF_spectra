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

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
outputdir <- args[2]
nsign <- as.numeric(args[3])
setwd(outputdir)
save.image("../mut_pat.Rdata")

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


garb <- sapply(base_clone, function(x) { 
  #pdf(paste0(x, ".96_profile.pdf"))
  ggsave(file=paste0(x, ".96_profile.pdf"), plot=plot_96_profile(mut_mat[,grepl(x, colnames(mut_mat))]))
  #graphics.off()
} )

### NMF for signatures
nrun_estimate <- 50
nrun <- 150
seed <- 123456
pdf("nmf_num_sign.pdf")
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:5, method="brunet", nrun=nrun_estimate, seed=seed)
plot(estimate)
graphics.off()

nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)

names_sign <- paste0("sign_", base, "_", seq(1, nsign))
colnames(nmf_res$signatures) <- names_sign
rownames(nmf_res$contribution) <- names_sign

pdf("signatures.96_profile.pdf")
plot_96_profile(nmf_res$signatures, condensed = TRUE)
graphics.off()
#plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "relative")

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

pdf("signatures_contribution_barplot.pdf")
mypc(nmf_res$contribution, nmf_res$signature, mode = "relative")
graphics.off()

pdf("signature_contribution.pdf")
plot_contribution_heatmap(nmf_res$contribution,sig_order = names_sign,cluster_samples=F)
graphics.off()


garb <- sapply(sample_names, function(x) { 
  #pdf(paste0(x, ".signature_evaluation.pdf"))
  pos <- which(sample_names==x)
  ggsave(file=paste0(x, ".signature_evaluation.pdf"), plot=plot_compare_profiles(mut_mat[,pos],nmf_res$reconstructed[,pos],profile_names = c("Original", "Reconstructed"),condensed = TRUE))
  #graphics.off()
} )



### Cosmic signatures
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
## why NA?

hclust_cosmic = cluster_signatures(cancer_signatures, method = "average") # store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
pdf("cosine_sign_cosmic.pdf")
plot_cosine_heatmap(cos_sim_samples_signatures,col_order = cosmic_order,cluster_rows = TRUE)
graphics.off()

ff <- fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(ff$contribution) > 10)
# Plot contribution barplot
pdf("cosmic_contribution_topsign.pdf") ## not pretty
plot_contribution(ff$contribution[select,],cancer_signatures[,select],coord_flip = FALSE,mode = "absolute")
graphics.off()
pdf("cosmic_contribution_heatmap.pdf")
plot_contribution_heatmap(ff$contribution,cluster_samples = TRUE,method = "complete")
graphics.off()



garb <- sapply(sample_names, function(x) { 
  #pdf(paste0("cosmic_",x,".signature_evaluation.pdf"))
  pos <- which(sample_names==x)
  ggsave(file=paste0("cosmic_",x,".signature_evaluation.pdf"), plot=plot_compare_profiles(mut_mat[,pos], ff$reconstructed[,pos],profile_names = c("Original", "Reconstructed"), condensed = TRUE))
  #graphics.off()
})

us <- nmf_res$signatures
colnames(us) <- names_sign
hclust_cosmic_us = cluster_signatures(cbind(us, cancer_signatures), method = "average") # store signatures in new order
pdf("hclust_cosmic_us.pdf")
plot(hclust_cosmic_us)
graphics.off()

cos_sim_ori_rec <- cos_sim_matrix(mut_mat, ff$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)

ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
  geom_bar(stat="identity", fill = "skyblue4") +
  coord_cartesian(ylim=c(0.8, 1)) +
  coord_flip(ylim=c(0.8,1)) +
  ylab("Cosine similarity\n original VS reconstructed") + 
  xlab("") + xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +theme_bw()+
  theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +geom_hline(aes(yintercept=.95))
ggsave("cosine_samples_cosmic.pdf")

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
# exit(0)
# chromosomes <- seqnames(get(ref_genome))[1:22]
# gb <- lapply(sample_names,function(x) {
#   ggsave(file=paste0(x,"_waterfall.png"), plot=plot_rainfall(vcf = vcfs[[x]], title=x, chromosomes=chromosomes))
#   })

