#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
infile_1 <- args[1]
infile_2 <- args[2]
name1 <- args[3]
name2 <- args[4]
column <- args[5]
prefix <- args[6]

data_1 <- read.table(gzfile(infile_1), sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)
data_2 <- read.table(gzfile(infile_2), sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)

n1 <- nrow(data_1)
n2 <- nrow(data_2)
clist <- intersect(data_1[,column], data_2[, column])
common <- length(clist)
cat(n1)
cat("\t")
cat(n2)
cat("\t")
cat(common)
cat("\n")


write.table(data.frame("ID"=clist), file = paste0(prefix, ".common.ids.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

library(VennDiagram)
venn.diagram(
  x = list(data_1[,column], data_2[,column]),
  category.names = c(name1, name2),
  filename = paste0(prefix, "_venn.png"),
  imagetype= "png",
  cat.cex = 0.8,
  cat.default.pos = "text"
)

