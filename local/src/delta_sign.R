library(MutationalPatterns)

load('/mnt/trcanmed/snaketree/prj/AF_spectra/dataset/mut_pat_signatures_vitro_2/mut_pat2.Rdata')
merged_contribution <- ff$contribution
load('/mnt/trcanmed/snaketree/prj/AF_spectra/dataset/MutationalPattern_bulk/mut_pat_signatures_2/mut_pat2.Rdata')
bulk_contribution <- ff$contribution
bulk_contribution <- t(ff$contribution)
merged_contribution <- t(merged_contribution)
merged_contribution_norm = merged_contribution/rowSums(merged_contribution)
bulk_contribution_norm = bulk_contribution/rowSums(bulk_contribution)

merged_contribution_norm <- merged_contribution_norm[match(rownames(bulk_contribution_norm), rownames(merged_contribution_norm)),]

contribution_norm <- merged_contribution_norm-bulk_contribution_norm

hc.sample = hclust(dist(contribution_norm), method = method)
method= "complete"
sample_order <- rownames(contribution_norm)[hc.sample$order]
sig_order = colnames(contribution_norm)
#https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt
contribution_norm.m = melt(contribution_norm)
colnames(contribution_norm.m) = c("Sample", "Signature",    "Contribution")
contribution_norm.m$Sample = factor(contribution_norm.m$Sample, levels = sample_order)
contribution_norm.m$Signature = factor(contribution_norm.m$Signature, levels = sig_order)
heatmap = ggplot(contribution_norm.m, aes(y = Signature, x = Sample, fill = Contribution, order = Sample)) + geom_tile(color = "white") + 
  scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "Relative contribution \n Difference", limits=c(-0.3,0.4) ) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15)) + labs(x = NULL, y = NULL)

#  heatmap = heatmap + geom_text(aes(label = round(Contribution, 
#                                                 2)), size = 3)
  dhc = as.dendrogram(hc.sample)
  ddata = dendro_data(dhc, type = "rectangle")
  dendrogram = ggplot(segment(ddata)) + geom_segment(aes(x = x,  y = y, xend = xend, yend = yend)) + scale_y_reverse(expand = c(0.2, 0)) + theme_dendro()
  cowplot::plot_grid(heatmap, dendrogram, align="v", ncol=1, rel_heights = c(1, 0.3), axis="lr")
  