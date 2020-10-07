library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg19'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)

files <- c("/home/egrassi/P1.merged.vcf.gz","/home/egrassi/P2.merged.vcf.gz","/home/egrassi/P3.merged.vcf.gz")
sample_names <- c('cl1',"cl2","cl3")


vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
plot_spectrum(type_occurrences, CT = TRUE, by=sample_names)


nrun_estimate <- 50
nrun <- 150
seed <- 123456
# we do not add 0.0001 to mut_mat because extract_signatures already do so (other fun used later TODO FIXME?)
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:4, method="brunet", nrun=nrun_estimate, seed=seed)
plot(estimate)
nsign <- 2 # from the plot optimal n. is 2 considering cophenetic...
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)


names_sign <- paste0("PDO-MA_", seq(1, nsign))
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

############### private of regions
path <- '/scratch/trcanmed/AF_spectra/local/share/data/clevers/private_clones/'
f <- c("P1T1.vcf.gz","P1T2.vcf.gz","P1T3.vcf.gz","P1T4.vcf.gz","P2T1.vcf.gz","P2T2.vcf.gz","P2T4.vcf.gz","P2T5.vcf.gz","P2T6.vcf.gz","P3T1.vcf.gz","P3T2.vcf.gz","P3T3.vcf.gz","P3T4.vcf.gz")

files <- paste0(path, f)
sample_names <- c("P1.T1.1","P1.T1.4","P1.T2.4","P1.T3.2","P1.T3.3","P1.T4.3","P1.T4.4","P2.T1.1","P2.T1.3","P2.T2.5","P2.T4.2","P2.T4.3","P2.T5.1","P2.T5.4","P2.T6.2","P2.T6.6","P3.T1.1","P3.T1.4","P3.T2.1","P3.T2.2","P3.T3.2","P3.T3.4","P3.T4.1","P3.T4.2")


vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)


nrun_estimate <- 50
nrun <- 150
seed <- 123456
# we do not add 0.0001 to mut_mat because extract_signatures already do so (other fun used later TODO FIXME?)
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:4, method="brunet", nrun=nrun_estimate, seed=seed)
plot(estimate)
nsign <- 3 
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)


names_sign <- paste0("PDO-MA_", seq(1, nsign))
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

## private clones

f <- c("P1.T1.1.vcf.gz","P1.T1.4.vcf.gz","P1.T2.4.vcf.gz","P1.T3.2.vcf.gz","P1.T3.3.vcf.gz","P1.T4.3.vcf.gz","P1.T4.4.vcf.gz","P2.T1.1.vcf.gz","P2.T1.3.vcf.gz","P2.T2.5.vcf.gz","P2.T4.2.vcf.gz","P2.T4.3.vcf.gz","P2.T5.1.vcf.gz","P2.T5.4.vcf.gz","P2.T6.2.vcf.gz","P2.T6.6.vcf.gz","P3.T1.1.vcf.gz","P3.T1.4.vcf.gz","P3.T2.1.vcf.gz","P3.T2.2.vcf.gz","P3.T3.2.vcf.gz","P3.T3.4.vcf.gz","P3.T4.1.vcf.gz","P3.T4.2.vcf.gz")
path <- '/scratch/trcanmed/AF_spectra/local/share/data/clevers/private_clones/'

files <- paste0(path, f)
sample_names <- c("P1.T1.1","P1.T1.4","P1.T2.4","P1.T3.2","P1.T3.3","P1.T4.3","P1.T4.4","P2.T1.1","P2.T1.3","P2.T2.5","P2.T4.2","P2.T4.3","P2.T5.1","P2.T5.4","P2.T6.2","P2.T6.6","P3.T1.1","P3.T1.4","P3.T2.1","P3.T2.2","P3.T3.2","P3.T3.4","P3.T4.1","P3.T4.2")


vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
plot_spectrum(type_occurrences, CT = TRUE, by=sample_names)


nrun_estimate <- 50
nrun <- 150
seed <- 123456
# we do not add 0.0001 to mut_mat because extract_signatures already do so (other fun used later TODO FIXME?)
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:4, method="brunet", nrun=nrun_estimate, seed=seed)
plot(estimate)
nsign <- 2 
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)


names_sign <- paste0("Clevers_singleclones_", seq(1, nsign))
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

