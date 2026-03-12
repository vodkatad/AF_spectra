# Libraries
if (!require(vcfR)) stop("Package 'vcfR' missing\n.")
if (!require(phangorn)) stop("Package 'phangorn' missing\n.")
if (!require(phytools)) stop("Package 'phytools' missing\n.")

#model <- snakemake@params[['model']]
# Where is the vaf?
tab_bin_path = snakemake@input[["tab_bin"]]
#vaf_path = '/home/data/work/stash/MA/CRC1307/platy/platypus_filtered.vcf.gz'
# What's the normal sample?
#normal  = snakemake@params[["normal"]]
# pdf output
pdf = snakemake@output[["pdf"]]
#pdf = "CRC1307_tree.pdf"
# nexus output
nexus.file = snakemake@output[["nexus"]]
#nexus.file ='CRC1307_nexus'
# Which patient?
patient = snakemake@params[["model"]] # keep to use for methods?
#patient = 'CRC1307'
# vaf threshold
#vt = snakemake@params[["vaf_threshold"]]
outimage = snakemake@output[['rimage']]

vafs_bin = read_table(tab_bin_path, sep="\t", header=T) # tsv with muts on rows and samples on columns, eg:
#              CRC1307-02-0 CRC1307-02-1-A CRC1307-02-1-B CRC1307-02-1-E
#chr1_10011343             0              0              0              0
#chr1_100189664            0              0              0              0
#chr1_100225993            1              1              1              1
#chr1_100247652            1              1              1              1
#chr1_100320358            0              0              0              0
#chr1_100425866            0              0              0              0

# Rownames
rownames(vafs_bin) = rownames(vafs)

# remove vivo if needed

#MA: note that I am not removing all 0 muts but an empiric check for 1502 showed that there are no differences (as expected) doing so, just slightly different bs percentages.
#ss <- rowSums(vafs_bin)
#length(ss)
#nss <- names(ss[ss==0])
#length(nss)
#v <- vafs_bin[!rownames(vafs_bin) %in% nss,]


# Make phyData
phyHs = as.phyDat(t(vafs_bin), type = "USER", levels = c("0", "1"))

# Bootstrap trees
bs = bootstrap.phyDat(phyHs, pratchet, bs = bs)
# Make the tree we'll go with
terry = pratchet(phyHs, maxit = 1000)
# Root it
#terry = root(terry, outgroup= normal, resolve.root = TRUE)
terry = root(terry, resolve.root = FALSE) # not checked what happens removing outgroup=normal and setting resolve.root=FALSE
# Get branch lengths
terry = acctran(terry, phyHs)
#Plot it with the bootstrap values
plot = plotBS(tree = terry, BStrees = bs, type = "phylogram", p = 0)

#Write out the tree as a nexus file
writeNexus(plot, file=nexus.file)
#Dan style trees
NexusTree = read.nexus(nexus.file)
if (class(NexusTree)=="multiPhylo"){
  NexusTree = NexusTree[[1]]
}

pdf(pdf)

#plot the tree

plotBS(NexusTree,
       BStrees = bs,
       edge.width=4,
       label.offset=4.0,
       font=0.5, 
       cex=1, 
       main = patient, 
       type = "phylogram")
add.scale.bar(lwd=8, cex=1.5)
dev.off()

save.image(outimage)
