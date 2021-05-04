load('/scratch/trcanmed/AF_spectra/dataset/CRC1502_clones_all/tree/tree_bulk_vitro.Rdata')
plotBS(NexusTree,
BStrees = bs,
edge.color=clcolr,
edge.width=3,
label.offset=4.0,
font=0.5,
cex=1,
main = patient,
type = "phylogram")


og  =  match(normal, NexusTree$tip.label)
oge =  NexusTree$edge[NexusTree$edge[,2]==og][1]

#BG colr vector
clcolr =  rep("#DF8476", dim(NexusTree$edge)[1])
#Change the internal edge colours
#clcolr[NexusTree$edge[,2] %in% NexusTree$edge[,1]] = '#F9DD4B'
#Add og colr
clcolr[NexusTree$edge[,1]==oge] = '#4999D3'

layout(matrix(1:3, 1, 3), widths=c(0.1,0.8,0.1))
plot.new()
plot(NexusTree, type="phylogram", edge.color=clcolr,
     edge.width=3,
     label.offset=4.0,
     font=1,
     cex=1.5,
     main = patient, font.main=2, cex.main=2)
#locator() 
add.scale.bar(x=47898.85, y= 1.1833010,lwd=2, cex=1.5)