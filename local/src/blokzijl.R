library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg19'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)

files <- c("C30913DSubclone109_C30913DBiopsy_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf","C30913Dsubclone125_C30913DBiopsy_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf","STE0072SC12A_BOXTELBLOOD0072_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf","STE0072SC23B_BOXTELBLOOD0072_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf","STE0072SC31A_BOXTELBLOOD0072_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf","STE0076SC12C_BOXTELBLOOD0076_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf" ,"STE0076SC23A_BOXTELBLOOD0076_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf","STE0076SC32A_BOXTELBLOOD0076_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf","subclone105_BIOPSY17513D_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf","subclone33_BIOPSY17513D_Q100_PASS_20X_autosomal_noSNP_nonRECUR_final_true_invitro.vcf")
files <- paste0('/scratch/trcanmed/AF_spectra/local/share/data/Blokzijl/', files)
sample_names <- seq(1, length(files))
sample_names <- c("Liver_14-d","Liver_14-e","Intestine_5-d","Intestine_5-e","Intestine_5-f","Intestine_6-d","Intestine_6-e","Intestine_6-f","Liver_17-d","Liver_17-e")
# egrassi@godot:/scratch/trcanmed/AF_spectra/local/share/data/Blokzijl$ tr -s " " "\t"<  meta  | awk '{print $3"_"$5}' | tr "\n" "," |sed 's/,/","/g'

vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
#plot_spectrum(type_occurrences, CT = TRUE, by=sample_names)


nrun_estimate <- 50
nrun <- 150
seed <- 123456
# we do not add 0.0001 to mut_mat because extract_signatures already do so (other fun used later TODO FIXME?)
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:4, method="brunet", nrun=nrun_estimate, seed=seed)
plot(estimate)
nsign <- 3 # from the plot optimal n. is 3 considering cophenetic...
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)


names_sign <- paste0("Blokzjil_", seq(1, nsign))
colnames(nmf_res$signatures) <- names_sign
rownames(nmf_res$contribution) <- names_sign

#compare with cosmic
cosmic <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "") # ???
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


us <- nmf_res$signatures
colnames(us) <- names_sign
hclust_cosmic_us = cluster_signatures(cbind(us, cancer_signatures), method = "average") 
plot(hclust_cosmic_us)

plot_contribution_heatmap(nmf_res$contribution)

### vs signal
signal <- '/home/egrassi/signal_colorectal.txt'
crc_signatures = read.table(signal, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), crc_signatures$Substitution.Type)
# Reorder cancer signatures dataframe>
crc_signatures = crc_signatures[as.vector(new_order),]
crc_signatures = as.matrix(crc_signatures[,2:(ncol(crc_signatures)-1)])
# Add trinucletiode changes names as row.names>
row.names(crc_signatures) = crc_signatures$Substitution.Type
ff <- fit_to_signatures(mut_mat, crc_signatures)
ff$contribution
which(rowSums(ff$contribution) > 10)
select <- which(rowSums(ff$contribution) > 10)
plot_contribution_heatmap(ff$contribution,cluster_samples = TRUE,method = "complete")

cos_sim_ori_rec <- cos_sim_matrix(mut_mat, ff$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)

ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
  geom_bar(stat="identity", fill = "skyblue4") +
  coord_cartesian(ylim=c(0.7, 1)) +
  coord_flip(ylim=c(0.7,1)) +
  ylab("Cosine similarity\n original VS reconstructed") + 
  xlab("") + xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +theme_bw()+
  theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +geom_hline(aes(yintercept=.95))