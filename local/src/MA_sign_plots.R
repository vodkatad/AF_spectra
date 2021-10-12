load('/scratch/trcanmed/AF_spectra/dataset/4_all_mutpat2_denovo/mutpat.Rdata')
names <- c('MA_MSI','MA_8','MA_5', 'MA_18')
colnames(nmf_res$signatures) <- names
rownames(nmf_res$contribution) <- names
pal <- c("#999999", "#E69F00", "#56B4E9",'#79d279')
mypc(nmf_res$contribution, nmf_res$signature, mode = "relative", palette=pal)

plot_96_profile(nmf_res$signatures[,2, drop=F], condensed = TRUE)
plot_96_profile(cancer_signatures[,8, drop=F], condensed=TRUE)
cos_sim_matrix(us, cancer_signatures)
coss <- cos_sim_matrix(us, cancer_signatures)
apply(coss, 1, max)
apply(coss, 2, max)
