dist <- snakemake@input[['dist']]

data <- read.table(dist, sep="\t", stringsAsFactors=FALSE)
colnames(data) <- gsub('.', '-', colnames(data), fixed=T)
save.image('pippo.Rdata')
names <- rownames(data)
time <- sapply(names, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})
model <- sapply(names, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
clone_begin <- sapply(names, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
clone_end <- sapply(names, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
# all t1 with mes0, all t2 with t1
maiscrivereifor <- function(i) {
     # we skip bulk and matched normals, with a variability of names
     if (grepl('LMO', names[i]) || grepl('PRO', names[i]) || grepl('H', names[i]) || grepl('diploid', names[i])) {
          return(c('', '', ''))
     }
     if (time[i] == '1') {
          wanted <- paste0(model[i], "-", clone_begin[i], '-0')
          add <- c(names[i], wanted, data[names[i], wanted])
          return(unlist(add))
     } else if (time[i] == '2') {
          clone_b <- substr(clone_begin[i], 1, 2)
          clone_e <- substr(clone_begin[i], 3, 3)
          wanted <- paste0(model[i], "-", clone_b, '-1-', clone_e)
          add <- c(names[i], wanted, data[names[i], wanted])
          return(unlist(add))
     }
     return(c('', '', ''))
}

res <- sapply(seq(1, length(names)), maiscrivereifor)
#res <- Filter(function(x) x[] == 3, res)
res <- as.data.frame(t(res), stringsAsFactor=FALSE)
colnames(res) <- c('from', 'to', 'dist')
res <- res[res$dist != "",]
res$dist <- as.numeric(res$dist)
write.table(res, file=snakemake@output[['tips']], sep="\t", quote=FALSE, row.names=FALSE)
