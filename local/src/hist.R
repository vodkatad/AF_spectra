library(ggplot2)
selection <- snakemake@input[[1]]
#out <- snakemake@output[[1]]
outall <- snakemake@output[['all']]
outtable <- snakemake@output[['table']]
outselected <- snakemake@output[['selected']]
depths <- snakemake@params[[1]][["depths"]]
lower <- snakemake@params[[1]][["lower"]]
upper <- snakemake@params[[1]][["upper"]]
save.image("p.RData")
afs <- vector("list", length(depths))
names(afs) <- depths
total <- 0
for (i in seq(2, length(snakemake@input))) {
  afs[[i-1]] <- read.table(snakemake@input[[i]], header=FALSE, sep="\t")
  colnames(afs[[i-1]]) <- c("ID","AF")
  total <- nrow(afs[[i-1]]) + total
}

if (length(depths) != length(afs)) {
  stop("Bad llama")
}
wanted <- read.table(selection, header = FALSE, sep="\t")
colnames(wanted) <- 'ID'


for (i in seq(1, length(afs))) {
  afs[[i]]$depth <- depths[i]
}
all <- data.frame()
for (i in seq(1, length(afs))) {
  all <- rbind(all, afs[[i]])
}

colnames(all) <- c("ID", "AF", "depth")
# reorder depths
all$depth <- as.factor(all$depth)
if (length(levels(all$depth)) == 7) {
    all$depth <- factor(all$depth, levels=c("10X","30X","50X","80X","100X","200X","300X"))
} else if (length(levels(all$depth)) == 8) {
    all$depth <- factor(all$depth, levels=c("10X","30X","50X","80X","100X","200X","300X","gs"))
} else {
    stop("unknown depths")
}

ggplot(all, aes(x=AF, color=depth))+geom_freqpoly(bins=50)
ggsave(outall)

selected <- all[all$ID %in% wanted$ID,]

ggplot(selected, aes(x=AF, color=depth))+geom_freqpoly(bins=50)
ggsave(outselected)

df <- as.data.frame(table(selected$depth))

sums <- data.frame()
#summaries
for (i in seq(1, length(depths))) {
  x <- selected[selected$depth == depths[i],]$AF
  seen <- length(x[x<=upper&x>=lower])
  sums <- rbind(sums, c(summary(x), sd(x), seen), deparse.level=0)
}
colnames(sums) <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.", "sd","seen")
rownames(sums) <- depths
rownames(df) <- df$Var1
df$Var1 <- NULL
save.image("p.Rdata")
m <- merge(df, sums, by="row.names")
rownames(m) <- m[,1]
m[,1] <- NULL
colnames(m)[1] <- 'nmuts'
n <- max(m$nmuts)
m$perc_found <- (m$nmuts/n)*100
m$perc_interval <- (m$seen/n)*100
m <- m[levels(all$depth),]

write.table(m, outtable, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
