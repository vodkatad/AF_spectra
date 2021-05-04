# Libraries
if (!require(vcfR)) stop("Package 'vcfR' missing\n.")
if (!require(phangorn)) stop("Package 'phangorn' missing\n.")
if (!require(phytools)) stop("Package 'phytools' missing\n.")

#model <- snakemake@params[['model']]
# Where is the vaf?
vaf_path = snakemake@input[["vcf"]]
bulk_vaf_path = snakemake@input[["bulk"]]
#vaf_path = '/home/data/work/stash/MA/CRC1307/platy/platypus_filtered.vcf.gz'
# What's the normal sample?
#normal  = snakemake@params[["normal"]]
normal = snakemake@params[['normal']]
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
outimage = snakemake@output[['rimage']]
vt = 0.1
bs = 10000
#bs = 1
vivi = snakemake@wildcards[['vv']]

chrs = snakemake@params[['chrs']]
keep_chr <- unlist(strsplit(chrs, ','))

save.image('pippo.Rdata')
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

# remove vivo if needed
if (vivi == "vitro") {
  cnames <- colnames(vafs_bin)
  vafs_bin <- vafs_bin[, !grepl('-M', cnames, fixed=T)]
}

# remove unwanted chrs
rn <- rownames(vafs_bin)
seen_chr <- sapply(strsplit(rn,"_"), function(x){x[[1]]})
vafs_bin <- vafs_bin[seen_chr %in% keep_chr,]

## same procedure for mutect vcf for bulk
# Read in the vcf
bvcf = read.vcfR(bulk_vaf_path)
# Calculate the vaf
bvafs = extract.gt(bvcf, element = "AF", as.numeric = TRUE)
# Use threshold for vaf
bvafs = bvafs >= vt
# Convert to binary table
bvafs_bin = apply(bvafs, 2, as.integer)
# Rownames
rownames(bvafs_bin) = rownames(bvafs)

# remove vivo if needed
if (vivi == "vitro") {
  cnames <- colnames(bvafs_bin)
  bvafs_bin <- bvafs_bin[, !grepl('-M', cnames, fixed=T)]
}


# remove unwanted chrs
rn <- rownames(bvafs_bin)
seen_chr <- sapply(strsplit(rn,"_"), function(x){x[[1]]})
bvafs_bin <- bvafs_bin[seen_chr %in% keep_chr,]


vafs_bin <- merge(vafs_bin, bvafs_bin, all=T, by=0)
vafs_bin[is.na(vafs_bin)] <- 0
rownames(vafs_bin) <- vafs_bin$Row.names
vafs_bin$Row.names <- NULL
save.image('pippo.Rdata')

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
       font=0.5, 
       cex=1, 
       main = patient, 
       type = "phylogram")
add.scale.bar(lwd=8, cex=1.5)
dev.off()

save.image(outimage)