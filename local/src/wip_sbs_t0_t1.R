load('/scratch/trcanmed/AF_spectra/datasetV2/p.Rdata')


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
  geom_hline(aes(yintercept=.95))+theme(legend.position="none")



data <- as.matrix(ff$contribution)
#select <- which(rowSums(data) > 10)
data = t(data)
data = data/rowSums(data)
#select <- which(colSums(data)>1)
#data <- data[,select]
annot_rows <- data.frame(row.names=rownames(data))

#annot_rows[grepl('M',rownames(annot_rows)), 'sample'] <- 'vivo'
annot_rows$time <-  as.factor(unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[[3]] })))
annot_rows$model <- unlist(lapply(strsplit(rownames(annot_rows),'-'), function(x){ x[1] }))
annot_rows$model <- as.factor(unlist(lapply(strsplit(annot_rows$model,'_'), function(x){ x[1] })))
#annot_cols <- data.frame(signature=c('Clock','Clock','ROS','smoke','MMR','MMR','NA','NA'), row.names=colnames(data))

# cosmic specific XXX WARNING
colnames(data) <- unlist(lapply(strsplit(colnames(data),'.', fixed=TRUE), function(x){ x[length(x)] }))

pheatmap(data, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE,  annotation_row=annot_rows, color=brewer.pal(9, 'PuBu'))

#
data2 <- data[!grepl('-M', rownames(data)),]
annot_rows2 <- annot_rows[annot_rows$time %in% c(0, 1),]
pheatmap(data2, fontsize_row = 9, fontsize_col=9, show_colnames = TRUE,  cluster_cols=FALSE,  annotation_row=annot_rows2, color=brewer.pal(9, 'PuBu'))
