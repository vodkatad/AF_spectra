setwd("/mnt/trcanmed/snaketree/prj/AF_spectra/dataset/CRC1307_platypus_nobin/")
library(MutationalPatterns)
load('/mnt/trcanmed/snaketree/prj/AF_spectra/local/share/data/mutational_patterns/data/MutPat_object.rds'
)
load('/mnt/trcanmed/snaketree/prj/AF_spectra/local/share/data/mutational_patterns/data/MutPat_object.rds'
)
readRDS('/mnt/trcanmed/snaketree/prj/AF_spectra/local/share/data/mutational_patterns/data/MutPat_object.rds')
MutPat_object = readRDS(paste(dir, "/data/MutPat_object.rds",sep=""))
#head(multipleClasses())
MutPat_object <- readRDS('/mnt/trcanmed/snaketree/prj/AF_spectra/local/share/data/mutational_patterns/data/MutPat_object.rds')
MutPat_object$vcf
MutPat_object$vcf$`13-b`
nrow(MutPat_object$vcf$`13-b`)
length(MutPat_object$vcf$`13-b`)
lapply(MutPat_object$vcf, length)
summary(lapply(MutPat_object$vcf, length))
summary(sapply(MutPat_object$vcf, length))
library('BSgenome.Hsapiens.UCSC.hg38', character.only = TRUE)
files <- list.files(path="./", pattern = "*gain.vcf.gz", recursive = F, full.names = F)
files
files <- list.files(path="./", pattern = "*gain.vcf.gz", recursive = F, full.names = F)
sample_names <- unlist(lapply(files, function(x) {strstplit(x,"_")[[1]][1]} ))
sample_names <- unlist(lapply(files, function(x) {strsplit(x,"_")[[1]][1]} ))
sample_names
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
warnings()
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
type_occurrences
plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)
plot_spectrum(type_occurrences, CT = TRUE)
plot_spectrum(type_occurrences, CT = TRUE, by=sample_names)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
head(mut_mat)
plot_96_profile(mut_mat)
plot_96_profile(mut_mat[,1])
plot_96_profile(mut_mat[,c(1,2)])
plot_96_profile(mut_mat[,c(1,2,3)])
plot_96_profile(mut_mat[,seq(4,12)])
head(mut_mat[,seq(4,12)]
)
mut_mat <- mut_mat + 0.0001
library(NMF)
estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)
plot(estimate)
nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 30)
head(nmf_res$signatures)
nmf_res <- extract_signatures(mut_mat, rank = 3, nrun = 30)
colnames(nmf_res$signatures) <- c("Signature 1", "Signature 2", "Signature 3")
rownames(nmf_res$contribution) <- c("Signature 1", "Signature 2", "Signature 3")
plot_96_profile(nmf_res$signatures, condensed = TRUE)
nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 100)
rownames(nmf_res$contribution) <- c("Signature 1", "Signature 2")
colnames(nmf_res$signatures) <- c("Signature 1", "Signature 2")
plot_96_profile(nmf_res$signatures, condensed = TRUE)
plot_contribution(nmf_res$contribution, nmf_res$signature,+                          mode = "relative")
plot_contribution(nmf_res$contribution, nmf_res$signature,                          mode = "relative")
plot_contribution(nmf_res$contribution, nmf_res$signature,                          mode = "absolute")

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
plot_contribution_heatmap(nmf_res$contribution,sig_order = c("Signature 1", "Signature 2"))
plot_contribution_heatmap(nmf_res$contribution,sig_order = c("Signature 1", "Signature 2"),cluster_samples=F)

plot_compare_profiles(mut_mat[,1],nmf_res$reconstructed[,1],profile_names = c("Original", "Reconstructed"),condensed = TRUE)

###
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
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


plot_cosine_heatmap(cos_sim_samples_signatures,col_order = cosmic_order,cluster_rows = TRUE)

ff <- fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(ff$contribution) > 10)
# Plot contribution barplot
plot_contribution(ff$contribution[select,],cancer_signatures[,select],coord_flip = FALSE,mode = "absolute")
plot_contribution_heatmap(ff$contribution,cluster_samples = TRUE,method = "complete")

plot_compare_profiles(mut_mat[,1], ff$reconstructed[,1],profile_names = c("Original", "Reconstructed"),condensed = TRUE)

us <- nmf_res$signatures
colnames(us) <- c("Our Sign 1", "Our Sign 2")
hclust_cosmic_us = cluster_signatures(cbind(us, cancer_signatures), method = "average") # store signatures in new order
plot(hclust_cosmic_us)

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
                                 

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes)

strand_counts <- strand_occurrences(mut_mat_s, by=sample_names)
strand_bias <- strand_bias_test(strand_counts)
plot_strand(strand_counts[grepl('CRC1307-02',strand_counts$group),], mode = "relative")

clones <- unique(unlist(lapply(sample_names, function(x) {xx <-strsplit(x,"-"); paste0(xx[[1]][1],'-',xx[[1]][2])} )))
for (c in clones) {
  p0 <- plot_strand(strand_counts[grepl(c,strand_counts$group),], mode = "relative")
  p1 <- plot_strand_bias(strand_bias[grepl(c,strand_bias$group),])
  do.call(grid.arrange, list(p0,p1))
}


####
chromosomes <- seqnames(get(ref_genome))[1:22] # Make a rainfall plot
plots <- lapply(names(vcfs), function(x) { plot_rainfall(vcfs[[x]], title = x, chromosomes = chromosomes, cex = 1.5, ylim = 1e+09) })
do.call(grid.arrange, plots)

###
#library(bedr)
#regions = bed_to_granges("../MutationalPatterns/Homo_sapiens.GRCh37.75_autosomal_exon_merged_sorted.bed")
regions <- import("../MutationalPatterns/Homo_sapiens.GRCh37.75_autosomal_exon_merged_sorted.bed")
regions_list <- GRangesList(regions)
names(regions_list) <- c('exons')
#probably import is ok?
surveyed <- import("../../local/share/data/CRC1307_clones_mutect/callable_covered.bed.gz")
surveyed_list <- rep(list(surveyed), length(vcfs))
# warnings on different chr we don't care
distr <- genomic_distribution(vcfs, surveyed_list, regions_list)

distr_testwhole <- enrichment_depletion_test(distr, by = rep("whole", length(vcfs)))
distr_test <- enrichment_depletion_test(distr, by=sample_names)

plot_enrichment_depletion(distr_test)
plot_enrichment_depletion(distr_testwhole)

regions2 <- import("../../local/share/data/CRC1307_clones_mutect/all_normal_chr.bed")
rg <- split(regions2, seqnames(regions2))
regions_list <- GRangesList(rg[[6]])

names(regions_list) <- c("chr6")

distr <- genomic_distribution(vcfs, surveyed_list, regions_list)

distr_testwhole <- enrichment_depletion_test(distr, by = rep("whole", length(vcfs)))
distr_test <- enrichment_depletion_test(distr, by=sample_names)

plot_enrichment_depletion(distr_test)
plot_enrichment_depletion(distr_testwhole)
# mh, not ok with residuals of chi sq?  maybe its not ok to work on chromosomes...

