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
library(RColorBrewer)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
outputheat <- args[2]
outputcosine <- args[3]
palette <- args[4]

# infile is tsv vfc / samplename_bulk|vivo|vitro
vcfs <- read.table(input, sep="\t", header=FALSE, stringsAsFactors = FALSE)
vcf_files <- vcfs$V1
sample_names <- vcfs$V2

vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)


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

## signal TODO  new version of package has it already inside?
## other cosmic version
###



ff <- fit_to_signatures(mut_mat, ref_signatures)
#ff$contribution
#which(rowSums(ff$contribution) > 10)
select <- which(rowSums(ff$contribution) > 10)
#plot_contribution_heatmap(ff$contribution,cluster_samples = TRUE,method = "complete")


cos_sim_ori_rec <- cos_sim_matrix(mut_mat, ff$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$s = row.names(cos_sim_ori_rec)
cos_sim_ori_rec$sample <- unlist(lapply(strsplit(cos_sim_ori_rec$s,"_"), function(x) x[length(x)]))
cos_sim_ori_rec$n <- seq(1, nrow(cos_sim_ori_rec))
ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=reorder(s, n), fill=sample)) +
  geom_bar(stat="identity") +
  coord_flip(ylim=c(0.7,1)) +
  ylab("Cosine similarity\n original VS reconstructed") + 
  xlab("") +theme_bw()+ theme(axis.text.y = element_text(size=10), axis.title.x=element_text(size=15))+
  geom_hline(aes(yintercept=.95)) +scale_fill_brewer(palette = "Dark2")

ggsave(outputcosine)

data <- as.matrix(ff$contribution)
#select <- which(rowSums(data) > 10)
data = t(data)
data = data/rowSums(data)
#select <- which(colSums(data)>1)
#data <- data[,select]
annot_rows <- data.frame(row.names=rownames(data))

#annot_rows[grepl('M',rownames(annot_rows)), 'sample'] <- 'vivo'
annot_rows$sample <-  as.factor(unlist(lapply(strsplit(rownames(annot_rows),'_'), function(x){ x[length(x)] })))
annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] }))
annot_rows$model <- as.factor(unlist(lapply(strsplit(annot_rows$model,'_'), function(x){ x[1] })))
#annot_cols <- data.frame(signature=c('Clock','Clock','ROS','smoke','MMR','MMR','NA','NA'), row.names=colnames(data))
col <- brewer.pal(3,'Dark2')
names(col) <- levels(annot_rows$sample)
# cosmic specific XXX WARNING
colnames(data) <- unlist(lapply(strsplit(colnames(data),'.', fixed=TRUE), function(x){ x[length(x)] }))

cbPalette2 <- unlist(strsplit(palette, ','))

names(cbPalette2) <- levels(annot_rows$model)
annot_colors <- list(sample=col, model=cbPalette2)
pheatmap(data, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE,  annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), filename=outputheat)

############# manual AIRC
save.image('sign.Rdata')

q()
load('sign.Rdata')
annot_rows <- annot_rows[order(annot_rows$sample, annot_rows$model),]
data2 <- data[match(rownames(annot_rows), rownames(data)), ]

pheatmap(data2, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE, annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), filename="sign_cosmic_order1.pdf")


annot_rows <- annot_rows[order(annot_rows$model, annot_rows$sample),]
data3 <- data[match(rownames(annot_rows), rownames(data)), ]
pheatmap(data3, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE, annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), filename="sign_cosmic_order2.pdf")

pheatmap(t(data3), fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), file="t_signature_cosmic_order2.pdf" )

#data4 <- data3[grepl('vitroMA', rownames(data3)),]
data4 <- data3[!grepl('2nd', rownames(data3)),]
#annotation_rows2 <- annot_rows
#annotation_rows2$sample <- NULL
pheatmap(data4, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE, annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), filename="signnozoom_cosmic_order2.pdf")
pheatmap(t(data4), fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE, annotation_col=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), file="t_signnozoom_cosmic_order2.pdf", width=11.7, height=8.3)



data4 <- data4[,colnames(data4) %in% c(1,6,8,18)]
pheatmap(data4, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE, cluster_rows=FALSE, annotation_row=annot_rows, annotation_colors = annot_colors,  color=brewer.pal(9, 'PuBu'), filename="signzoom_cosmic_order2.pdf")


#
d2 <- data2[rownames(data2)=='CRC1307_bulk', colnames(data2) %in% c(1,8)]
d1 <- data2[rownames(data2)=='CRC1078_bulk', colnames(data2) %in% c(1,8)]
d <- as.data.frame(rbind(d1, d2))
d$names <- c('M1', 'M2')
l <- melt(d)
ggplot(data=l, aes(x=names, y=value, fill=variable))+geom_col(position="dodge")+theme_bw()+scale_fill_manual(values=c('darkgoldenrod','darkgreen'))

ddd <- as.data.frame(data2[, '8']/data2[, '1'])
colnames(ddd) <- c('SBS8_SBS1')
ddd <- ddd[grepl('bulk',rownames(ddd)),, drop=F]
ddd$model <- sapply(strsplit(rownames(ddd),"_"), function(x){x[[1]]})
ggplot(data=ddd, aes(y=SBS8_SBS1, x=model))+geom_col(position="dodge")+theme_bw()

colors <- "#cc3300,#f607b9,#9900ff,#155d00,#77a003,#0829fc,#ff9900,#ffff00"
cbPalette <- unlist(strsplit(colors, ','))


ggplot(data=ddd, aes(y=SBS8_SBS1, x=model, fill=model))+geom_col(position="dodge")+theme_bw()+scale_fill_manual(values=cbPalette)+scale_y_log10()

colors <- "#155d00,#77a003,#ff9900,#ffff00"
cbPalette <- unlist(strsplit(colors, ','))

ddd2 <- ddd[ddd$model %in% c('CRC1078', 'CRC1307', 'CRC1599LM', 'CRC1599PR'),]
ggplot(data=ddd2, aes(y=SBS8_SBS1, x=model, fill=model))+geom_col(position="dodge")+theme_bw()+scale_fill_manual(values=cbPalette)+scale_y_log10()
