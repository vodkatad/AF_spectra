# Libraries
if (!require(vcfR)) stop("Package 'vcfR' missing\n.")
if (!require(phangorn)) stop("Package 'phangorn' missing\n.")
if (!require(phytools)) stop("Package 'phytools' missing\n.")

#model <- snakemake@params[['model']]
# Where is the vaf?
vaf_path = snakemake@input[["vcf"]]
#vaf_path = '/home/data/work/stash/MA/CRC1307/platy/platypus_filtered.vcf.gz'
# What's the normal sample?
#normal  = snakemake@params[["normal"]]
normal = 'normal'
# pdf output
pdf = snakemake@output[["pdf"]]
#pdf = "CRC1307_tree.pdf"
# nexus output
nexus.file = snakemake@output[["nexus"]]
#nexus.file ='CRC1307_nexus'
# Which patient?
patient = snakemake@params[["model"]]
#patient = 'CRC1307'
# vaf threshold
#vt = snakemake@params[["vaf_threshold"]]
outimage = snakemake@input[['rimage']]
vt = 0.1
#bs = 1000 #10000
bs = 10000

# Read in the vcf
vcf = read.vcfR(vaf_path)
# What is the depth?
depth     = extract.gt(vcf, element = "NR", as.numeric = TRUE)
# How many reads support the variant?
variant   = extract.gt(vcf, element = "NV", as.numeric = TRUE)
# Calculate the vaf
vafs = variant / depth
# Use threshold for vaf
vafs = vafs >= vt
# Convert to binary table
vafs_bin = apply(vafs, 2, as.integer)
# Rownames
rownames(vafs_bin) = rownames(vafs)
# add fake outgroup/nrormal
vafs_bin <- cbind(rep(0, nrow(vafs_bin)),vafs_bin)
colnames(vafs_bin)[1] <- normal

# Make phyData
phyHs = as.phyDat(t(vafs_bin), type = "USER", levels = c("0", "1"))

# Bootstrap trees
bs = bootstrap.phyDat(phyHs, pratchet, bs = bs)
# Make the tree we'll go with
terry = pratchet(phyHs, maxit = 1000)
# Root it
terry = root(terry, outgroup= normal, resolve.root = TRUE)
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

#Identify the outgroup edge, to colour
og  =  match(normal, NexusTree$tip.label)
oge =  NexusTree$edge[NexusTree$edge[,2]==og][1]

#BG colr vector
clcolr =  rep("#DF8476", dim(NexusTree$edge)[1])
#Change the internal edge colours
clcolr[NexusTree$edge[,2] %in% NexusTree$edge[,1]] = '#F9DD4B'
#Add og colr
clcolr[NexusTree$edge[,1]==oge] = '#4999D3'
pdf(pdf)

#plot the tree

plotBS(NexusTree,
       BStrees = bs,
       edge.color=clcolr,
       edge.width=4,
       label.offset=4.0,
       font=1, 
       cex=1.3, 
       main = patient, 
       type = "phylogram")
add.scale.bar(lwd=8, cex=1.5)
dev.off()

save.image(outimage)