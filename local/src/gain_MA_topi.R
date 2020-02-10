#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
len <- args[1]
calls <- args[2]
gens <- args[3]
save.image("p.RData")

dcalls <- read.table(gzfile(calls), header=FALSE, stringsAsFactors=FALSE, sep="\t")
dlen <- read.table(len, header=FALSE)
glen <- dlen$V1
generation <- read.table(gens, header=TRUE, sep="\t", colClasses=c("character", "numeric", "character", "numeric"))

data <- dcalls[,c("V1","V16"), drop=F]
colnames(data) <- c("id", "class")
df <- as.data.frame(table(data$class))
df[df$Var1=="gain","Freq"]-df[df$Var1=="loss","Freq"]
gained <- df[df$Var1=="gain","Freq"]-df[df$Var1=="loss","Freq"]
name_num <- strsplit(len,"_")[[1]][1]
name_den <- strsplit(len,"_")[[1]][2]
name_den <- strsplit(name_den,'.', fixed=TRUE)[[1]][1]
clinfo <- strsplit(name_num, '-')[[1]]
clone <- clinfo[2]

clone_gen <- paste0(clone,'-', clinfo[3])
colnames(generation) <- c("CAMPIONI", "GENS", "TIPO","GIORNI")
generation <- generation[generation$CAMPIONI==clone_gen,]

days <- unique(generation$GIORNI)
if (length(days) != 1) {
    stop("qualquadranoncosa")
}

get_MR <- function(ngen, nmut, len) {
    return(nmut / (ngen*len))
}

get_MR_time <- function(days, nmut, len) {
    return(nmut / (days*glen))
}

MR <- sapply(unique(generation$TIPO), function(x ) get_MR(generation[generation$TIPO==x,"GENS"], gained, glen))
MR_time <- get_MR_time(days, gained, glen)

type <- strsplit(name_num, '-')[[1]][3]

res <- data.frame(clone=clone, end=name_num, start=name_den, class=type, n_gained=gained, len_cnok=glen, MR_conte=MR[1], MR_EDU=MR[1], MR_giorni=MR_time)

write.table(res, args[4], sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
