library(phytools)
library(usedist)
library(adephylo)

model <- snakemake@params[['model']]
infile <- snakemake@input[['tree']]

if (snakemake@params[['nexus']] == "yes") {
     tree <- read.nexus(infile)
} else {
     tree <- read.tree(infile)
}
cna_dists <- distTips(tree, tips=tree$tip.label)
#attr(mut_dists, "Labels") = gsub("-", "_", attr(mut_dists, "Labels"))
#mut_dists_sub = dist_subset(mut_dists, idx = attr(cna_dists, "Labels")[attr(cna_dists, "Labels") %in% attr(mut_dists, "Labels")])

#df = data.frame(Mutation_Phylogenetic_Distances = as.vector(mut_dists_sub), CNA_Phylogenetic_Distances = as.vector(cna_dists_sub))
res <- as.data.frame(as.matrix(cna_dists))
colnames(res) = gsub("_", "-", colnames(res))
rownames(res) = gsub("_", "-", rownames(res))
if (snakemake@params[['nexus']] != "yes") {
     rownames(res) <- paste0(model, "-", rownames(res))
     colnames(res) <- paste0(model, "-", colnames(res))
}
write.table(res, file=snakemake@output[['dist']], sep="\t", quote=FALSE)
