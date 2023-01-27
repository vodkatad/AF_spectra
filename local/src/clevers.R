
conflictRules("DelayedArray", exclude = "seed")
library(NMF)
library(SummarizedExperiment)

library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg19'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
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
path <- '/scratch/trcanmed/AF_spectra/local/share/data/clevers/private_sections/'
f <- c("P1T1.vcf.gz","P1T2.vcf.gz","P1T3.vcf.gz","P1T4.vcf.gz","P2T1.vcf.gz","P2T2.vcf.gz","P2T4.vcf.gz","P2T5.vcf.gz","P2T6.vcf.gz","P3T1.vcf.gz","P3T2.vcf.gz","P3T3.vcf.gz","P3T4.vcf.gz")

files <- paste0(path, f)
#sample_names <- c("P1.T1.1","P1.T1.4","P1.T2.4","P1.T3.2","P1.T3.3","P1.T4.3","P1.T4.4","P2.T1.1","P2.T1.3","P2.T2.5","P2.T4.2","P2.T4.3","P2.T5.1","P2.T5.4","P2.T6.2","P2.T6.6","P3.T1.1","P3.T1.4","P3.T2.1","P3.T2.2","P3.T3.2","P3.T3.4","P3.T4.1","P3.T4.2")
sample_names <-  sapply(f, function(x) {y<-strsplit(x, '.')[[1]][1]; return(y[1])})
sample_names <- names(sample_names)

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

f <- c("P1.T1.1.vcf.gz","P1.T1.4.vcf.gz","P1.T2.4.vcf.gz","P1.T3.2.vcf.gz","P1.T3.3.vcf.gz","P1.T4.3.vcf.gz","P1.T4.4.vcf.gz","P2.T1.1.vcf.gz","P2.T1.3.vcf.gz","P2.T2.5.vcf.gz","P2.T4.2.vcf.gz","P2.T4.3.vcf.gz","P2.T5.1.vcf.gz","P2.T5.4.vcf.gz","P2.T6.2.vcf.gz","P2.T6.6.vcf.gz")#,"P3.T1.1.vcf.gz","P3.T1.4.vcf.gz","P3.T2.1.vcf.gz","P3.T2.2.vcf.gz","P3.T3.2.vcf.gz","P3.T3.4.vcf.gz","P3.T4.1.vcf.gz","P3.T4.2.vcf.gz")
path <- '/scratch/trcanmed/AF_spectra/local/share/data/clevers/private_clones/'

files <- paste0(path, f)
sample_names <- c("P1.T1.1","P1.T1.4","P1.T2.4","P1.T3.2","P1.T3.3","P1.T4.3","P1.T4.4","P2.T1.1","P2.T1.3","P2.T2.5","P2.T4.2","P2.T4.3","P2.T5.1","P2.T5.4","P2.T6.2","P2.T6.6")#,"P3.T1.1","P3.T1.4","P3.T2.1","P3.T2.2","P3.T3.2","P3.T3.4","P3.T4.1","P3.T4.2")


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


# dio cane https://support.bioconductor.org/p/110844/

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

## trees


f <- c("P2.T1.1.vcf","P2.T1.3.vcf","P2.T2.5.vcf","P2.T4.2.vcf","P2.T4.3.vcf","P2.T5.1.vcf","P2.T5.4.vcf","P2.T6.2.vcf","P2.T6.6.vcf","P3.T1.1.vcf","P3.T1.4.vcf","P3.T2.1.vcf","P3.T2.2.vcf","P3.T3.2.vcf","P3.T3.4.vcf","P3.T4.1.vcf","P3.T4.2.vcf",
       "P1.T1.1.vcf","P1.T1.4.vcf","P1.T2.4.vcf","P1.T3.2.vcf","P1.T3.3.vcf","P1.T4.3.vcf","P1.T4.4.vcf")

path <- '/scratch/trcanmed/AF_spectra/local/share/data/clevers/trees/'

files <- paste0(path, f)
sample_names <- c("P2.T1.1","P2.T1.3","P2.T2.5","P2.T4.2","P2.T4.3","P2.T5.1","P2.T5.4","P2.T6.2","P2.T6.6","P3.T1.1","P3.T1.4","P3.T2.1","P3.T2.2","P3.T3.2","P3.T3.4","P3.T4.1","P3.T4.2", "P1.T1.1","P1.T1.4","P1.T2.4","P1.T3.2","P1.T3.3","P1.T4.3","P1.T4.4")

vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
plot_spectrum(type_occurrences, CT = TRUE, by=sample_names)

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

ff <- fit_to_signatures(mut_mat, cancer_signatures)
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

### better heatmap

data <- as.matrix(ff$contribution)
data = t(data)
data = data/rowSums(data)
annot_rows <- data.frame(row.names=rownames(data))

annot_rows$patient <- unlist(lapply(strsplit(rownames(annot_rows),'.', fixed=T), function(x){ x[1] }))
#annot_rows$clone <- as.factor(unlist(lapply(strsplit(rownames(annot_rows),'.', fixed=T), function(x){ paste0(x[2], '.', x[3]) })))

#col <- brewer.pal(3,'Dark2')
#names(col) <- levels(annot_rows$sample)
# cosmic specific XXX WARNING
colnames(data) <- unlist(lapply(strsplit(colnames(data),'.', fixed=TRUE), function(x){ x[length(x)] }))

#cbPalette2 <- unlist(strsplit(palette, ','))

#names(cbPalette2) <- levels(annot_rows$model)
#annot_colors <- list(sample=col, model=cbPalette2)

pheatmap(t(data), fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE,  annotation_col=annot_rows,   color=brewer.pal(9, 'PuBu'), width=6, height=8.3)


pheatmap(t(data), fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE,  annotation_col=annot_rows,   color=brewer.pal(9, 'PuBu'), width=6, height=8.3, filename="clevers.pdf")


pc <- 0.0001
ratiodf <- data.frame(row.names=rownames(data), ratio = log10((data[,'8']+pc)/(data[,'1']+pc)))

ratiodf$mut <- ifelse(unlist(lapply(strsplit(rownames(ratiodf),'_', fixed=T), function(x){ x[1] }))!="common", 'Private', 'Shared')
ratiodf$model <- unlist(lapply(strsplit(rownames(ratiodf),'_', fixed=T), function(x){ x[2] }))
ratiodf$id <- rownames(ratiodf)



ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

ggplot(data=ratiodf, aes(x=reorder(id, as.numeric(as.factor(model))), y=ratio, fill=mut))+geom_col(position="dodge")+
  theme_bw()+ctheme+
  xlab('Sample')+ylab("Log10 ratio") + ggtitle("Signature 8 / Signature 1")+
  scale_fill_manual(values=c('darkgoldenrod', 'darkgreen'))

ggsave('clevers_ratio.pdf', width=8, height=5, units="in")

