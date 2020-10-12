library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)
library(pheatmap)

# getting list of vcfs
f <- '/scratch/trcanmed/AF_spectra/dataset/AF_spectra/gained_clones_vitro'

f_comma <- readLines(f)
files <- unlist(strsplit(f_comma, ','))

# assembling samples names
pos <- regexpr("CRC\\d+\\-", files, perl=TRUE)
sample_names <- mapply(function(x, p) {substr(x,p,p+13)}, files, pos)

#ls CRC1307_platypus_nobin/*gain.vcf.gz CRC0282/platypus_nobin/*gain.vcf.gz | grep M | tr "\n" ","  > vitro_clones_gained

fv <- '/scratch/trcanmed/AF_spectra/dataset/vitro_clones_gained'
fv_comma <- readLines(fv)

filesv <- unlist(strsplit(fv_comma, ','))
filesv <- paste0('/scratch/trcanmed/AF_spectra/dataset/', filesv)
# assembling samples names
pos <- regexpr("CRC\\d+\\-", filesv, perl=TRUE)
sample_namesv <- mapply(function(x, p) {substr(x,p,p+14)}, filesv, pos)


files <- c(files, filesv)
sample_names <- c(sample_names, sample_namesv)

# starting with mut pat
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
# from tcga multiple alternative alleles/indels are removed, some samples (MSI?) many more indels (15-30 vs 1000+)

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

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
cos_sim_ori_rec$s = row.names(cos_sim_ori_rec)
cos_sim_ori_rec$sample <- 'vitro'
cos_sim_ori_rec[grepl('M', rownames(cos_sim_ori_rec)),'sample'] <- 'vivo'
colo <- brewer.pal(3, 'Dark2')[-1]
ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=s, fill=sample)) +
  geom_bar(stat="identity") +
  coord_flip(ylim=c(0.7,1)) +
  ylab("Cosine similarity\n original VS reconstructed") + 
  xlab("") +theme_bw()+ theme(axis.text.y = element_text(size=10), axis.title.x=element_text(size=15))+
  geom_hline(aes(yintercept=.95)) +scale_fill_manual(values=colo)
ggsave('bulk_clones_cosinesim.svg')

data <- as.matrix(ff$contribution)
data = t(data)
data = data/rowSums(data)
annot_rows <- data.frame(row.names=rownames(data))
#annot_rows[grepl('M',rownames(annot_rows)), 'sample'] <- 'vivo'
annot_rows$sample <- as.factor(ifelse(grepl('M',rownames(annot_rows)),'vivo', 'vitro'))
annot_rows$model <- as.factor(unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] })))
#annot_cols <- data.frame(signature=c('Clock','Clock','ROS','smoke','MMR','MMR','NA','NA'), row.names=colnames(data))
col <- brewer.pal(3,'Dark2')[-1]
names(col) <- levels(annot_rows$sample)

cbPalette2 <- c("#ff5733", 
                #"#9d01fc"
                "#f607b9",
                "#155d00",
                "#77a003",
                "#0829fc"
)
names(cbPalette2) <- levels(annot_rows$model)
annot_colors <- list(sample=col, model=cbPalette2)
pheatmap(data, fontsize_row = 9, fontsize_col=13, show_colnames = TRUE,  cluster_cols=FALSE, angle_col = 315,  annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(50, 'PuBu'))


# NMF de novo signature extraction
nrun_estimate <- 50
nrun <- 150
seed <- 123456
# we do not add 0.0001 to mut_mat because extract_signatures already do so (other fun used later TODO FIXME?)
mut_mat2 <- mut_mat + 0.0001
estimate <- nmf(mut_mat2, rank=2:5, method="brunet", nrun=nrun_estimate, seed=seed)

nsign <- 3 # from the plot optimal n. is 2 considering cophenetic...
nmf_res <- extract_signatures(mut_mat, rank = nsign, nrun = nrun)
# go on with 4 and evaluate similarities..

names_sign <- paste0("MA_", seq(1, nsign))
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

hclust_signal_us = cluster_signatures(cbind(us, crc_signatures), method = "average") 
plot(hclust_signal_us)

ff <- fit_to_signatures(mut_mat,  nmf_res$signatures)
ff$contribution
plot_contribution_heatmap(ff$contribution,cluster_samples = TRUE,method = "complete")


signal <- '/home/egrassi/signal_refsig.txt'
acrc_signatures = read.table(signal, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), acrc_signatures$Substitution.Type)
# Reorder cancer signatures dataframe> 
acrc_signatures = acrc_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names>
row.names(acrc_signatures) = acrc_signatures$Substitution.Type
# Keep only 96 contributions of the signatures in matrix
acrc_signatures = as.matrix(acrc_signatures[,2:ncol(acrc_signatures)])
hclust_asignal_us = cluster_signatures(cbind(us, acrc_signatures), method = "average") 
plot(hclust_asignal_us)
