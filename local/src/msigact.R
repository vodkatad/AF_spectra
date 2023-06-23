#!/opt/R/R-4.1.1/bin/Rscript
library(ggplot2)
library(ICAMS)
library(mSigAct)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
outputheat <- args[2]
outputcosine <- args[3]
palette <- args[4]
log_f <- args[5]

# infile is tsv vfc / samplename_bulk|vivo|vitro
vcfs <- read.table(input, sep="\t", header=FALSE, stringsAsFactors = FALSE)
vcf_files <- vcfs$V1
sample_names <- vcfs$V2

# I wanted to use expand in the rule that generates signinput_..., therefore single samples
# will create empty vcf if they do not have any in vitro/in vivo clone and we need to get rid of them here
keep_noempty <- file.size(vcf_files) != 0L
vcf_files <- vcf_files[keep_noempty]
sample_names <- sample_names[keep_noempty]

save.image('pippo.Rdata')
sink(log_f)
print("N samples:")
print(length(vcf_files))
table((unlist(lapply(strsplit(sample_names,'_'), function(x){ x[length(x)] }))))
print("2nd round")
print(table(grepl('2nd_', sample_names)))
print(vcf_files)
sink()

catalogs <- VCFsToCatalogs(vcf_files, ref.genome = "hg38",
                           variant.caller = "unknown", region = "genome", filter.status="PASS", names.of.VCFs	= sample_names)


spectra_from_cat <- catalogs$catSBS96[,seq(1, length(vcf_files)) , drop=FALSE]
attr(spectra_from_cat, 'class') <- c('matrix', 'array')

sbs8_presence <- SignaturePresenceTest(spectra = spectra_from_cat,
                                       sigs = cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96,
                                       target.sig.index = "SBS8",
                                       seed = 42,
                                       mc.cores = 4)

sbs1_presence <- SignaturePresenceTest(spectra = spectra_from_cat,
                                       sigs = cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96,
                                       target.sig.index = "SBS1",
                                       seed = 42,
                                       mc.cores = 4)


pd <- as.data.frame(sapply(sbs8_presence, function(x) {x$chisq.p}))
pd2 <- as.data.frame(sapply(sbs1_presence, function(x) {x$chisq.p}))
pd$sbs <- 'SBS8'
pd2$sbs <- 'SBS1'
colnames(pd) <- 'pval'
colnames(pd2) <- 'pval'
ppd <- rbind(pd, pd2)
colnames(ppd) <- c('pval', 'sbs')
ppd$sample <- c(rownames(pd), rownames(pd2))

 function(x){ x[1] })
write.table(ppd, file=outputcosine, sep="\t", quote=FALSE)

ppd$kind <-  as.factor(unlist(lapply(strsplit(ppd$sample,'_'), function(x){ x[length(x)] })))
ppd$model <- unlist(lapply(strsplit(ppd$sample,'_'), function(x){ x[1] }))

ppd$adjpval <- p.adjust(ppd$pval)
p <- function(x, data) {
  p1 <- data[data$kind==x,]
  return(ggplot(data=p1)+geom_col(aes(x=sample,y=-log10(adjpval), fill=sbs), position="dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(x))
}


plots <- lapply(levels(ppd$kind), p, ppd)
ggarrange(plotlist = plots, ncol=3)
ggsave(outputheat)
