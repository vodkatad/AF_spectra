library(rtracklayer)
#library(corrplot)
setwd("/mnt/trcanmed/snaketree/prj/AF_spectra/dataset/CRC1307_platypus_nobin/")
d <- read.table('all_gained.tsv', sep="\t", header=F, stringsAsFactors = F)
l <- unlist(strsplit(d$V1,':'))
dd <- as.data.frame(do.call(rbind, l))
colnames(dd) <- c('chr','start','ref','alt')

regions <- import('union_callable.bed.gz', format="BED")

lens <- data.frame(len=width(regions), names=seqnames(regions))
totlen <- sapply(levels(lens$names),  function(x) { l <- lens[lens$names==x,'len']; return(sum(l)) } )

# but we need to consider cn? :(
totmut <- nrow(dd)
chrs <- names(totlen)
total <- sum(totlen)
expected <- sapply(chrs, function(x) { t <- totlen[x]; names(t) <- NULL; (totmut / total)*t })
observed <- sapply(chrs, function(x) { f <- dd[dd$chr==x,]; nrow(f) })

m <- merge(as.data.frame(expected), as.data.frame(observed), by="row.names")

chisq <- sum( (m$o -m$e)**2 / m$e )
p <- pchisq(chisq, nrow(m), lower.tail=FALSE)

#http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
std_residuals <- data.frame(chr =m$Row.names, residuals=(m$o -m$e) / sqrt(m$e))

#corrplot(std_residuals, is.cor = FALSE)
std_residuals$chr <- as.character(std_residuals$chr)
std_residuals$chr <- as.numeric(substr(std_residuals$chr, 4, 4+nchar(std_residuals$chr)))
std_residuals <- std_residuals[order(std_residuals$chr),]
ggplot(std_residuals, aes(x=chr,y=residuals))+geom_bar(stat="identity", fill="steelblue")+ scale_x_continuous("chr", labels = as.character(std_residuals$chr), breaks = std_residuals$chr)+theme_bw()+ggtitle(paste0('ChiSQ residuals, pvalue=',p))