
if (!require(phangorn)) stop("Package 'phangorn' missing\n.")
if (!require(phytools)) stop("Package 'phytools' missing\n.")

image_path <- snakemake@input[["rimage"]]
myplot <- snakemake@output[["pdf"]]
load(image_path)
#load('/scratch/trcanmed/AF_spectra/dataset/CRC1502_clones_all/tree/tree_bulk_vitro.Rdata')

og  =  match(normal, NexusTree$tip.label)
oge =  NexusTree$edge[NexusTree$edge[,2]==og][1]

NexusTree$tip.label <- gsub(paste0(patient, '-'), '', NexusTree$tip.label)
#BG colr vector
clcolr =  rep("#DF8476", dim(NexusTree$edge)[1])
#Change the internal edge colours
#clcolr[NexusTree$edge[,2] %in% NexusTree$edge[,1]] = '#F9DD4B'
#Add og colr
clcolr[NexusTree$edge[,1]==oge] = '#4999D3'

pdf(myplot, width=3.5, height=2.7, pointsize=5)
#layout(matrix(1:3, 1, 3), widths=c(0.1,0.8,0.1))
#par(mar=c(1,1,1,1))
#plot.new()
plot(NexusTree, type="phylogram", edge.color=clcolr,
     edge.width=1,
     label.offset=0,
     font=1,
     cex=1,
     main = patient, font.main=2, cex.main=1)
#locator() 
#add.scale.bar(x=47898.85, y= 1.1833010,lwd=2, cex=1.5)
add.scale.bar(lwd=1, cex=1)
dev.off()