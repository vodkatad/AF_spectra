#!/usr/bin/env Rscript
library(MutationalPatterns)
library(gridExtra)
library(ggplot2)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
library(ref_genome, character.only = TRUE)

args <- commandArgs(trailingOnly = TRUE)
input_vcf <- args[1]
input_bed <- args[2]
input_bed_name <- args[3]
input_covered <- args[4]
input_cnv_dir  <- args[5] # put a warning in the snakefile about this input built inside the Rcode!
output_plot <- args[6]
output_tsv <- args[7]
save.image(paste0(output_tsv,".Rdata"))

files <- unlist(strsplit(input_vcf, ','))
whole_names <- unlist(lapply(files, function(x) {y <- basename(x);  strsplit(y,".", fixed=TRUE)[[1]][1] }))
sample_names <- unlist(lapply(files, function(x) {y <- basename(x); strsplit(y,"_")[[1]][1]} ))
vcfs <- read_vcfs_as_granges(files, sample_names, ref_genome)
base_clone <- unique(sapply(sample_names, function(x) {y<-strsplit(x, '-')[[1]]; return(paste0(y[1],'-',y[2]))}))
base <- unique(sapply(sample_names, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])}))

my_genomic_distribution = function (vcf_list, surveyed_list, region_list, cnv_list = NULL) 
{
  if (length(vcf_list) != length(surveyed_list)) 
    stop("vcf_list and surveyed_list must have the same length")
  if (is.null(names(region_list))) 
    stop(paste("Please set the names of region_list using:", 
               "    names(region_list) <- c(\"regionA\", \"regionB\", ...)", 
               sep = "\n"))
  df = data.frame()
  for (j in 1:length(region_list)) {
    for (i in 1:length(vcf_list)) {
      res = my_intersect_with_region(vcf_list[[i]], surveyed_list[[i]], 
                                  region_list[[j]], cnv_list[[i]])
      #> distr
      #region          sample n_muts surveyed_length         prob surveyed_region_length  expected observed
      #1   exons  CRC1307-02-1-A   2507      2818955411 8.893365e-07              114961572 102.23952      104
      #2   exons  CRC1307-02-1-B   2519      2818955411 8.935934e-07              114961572 102.72890      104
      # will need to change surveyed_length and surveyed_region_length according to CNVs
      res$region = names(region_list)[j]
      res$sample = names(vcf_list)[i]
      res = res[, c(7, 8, 1:6)]
      df = rbind(df, res)
    }
  }
  df$region = factor(df$region, levels = names(region_list))
  return(df)
}

my_intersect_with_region = function(vcf, surveyed, region, cnv=NULL) {
  # Number of mutations in vcf file
  n_muts = length(vcf)
  

  # Number of base pairs that were surveyed
  if (is.null(cnv)) {
    surveyed_length = sum(as.numeric(width(surveyed)))
  } else {
    # this is ok cause our cnvs are already overlapped with callable regions
    ov = findOverlaps(surveyed, cnv, ignore.strand = TRUE) # query, subject
    surveyed = cnv[subjectHits(ov)]
    surveyed_length = sum(as.numeric(width(surveyed)) * as.numeric(mcols(surveyed)$cn))
  }
  
  # Check if chromosome names are the same in the objects
  if (seqlevelsStyle(vcf) != seqlevelsStyle(surveyed))
    stop(paste("The chromosome names (seqlevels) of the VCF and the",
               "surveyed GRanges object do not match."))
  
  if (seqlevelsStyle(region) != seqlevelsStyle(surveyed))
    stop(paste("The chromosome names (seqlevels) of the surveyed and",
               "the region GRanges object do not match."))
  
  if (!is.null(cnv) && seqlevelsStyle(vcf) != seqlevelsStyle(cnv))
      stop(paste("The cnv names (seqlevels) of the VCF and the",
                 "CNV Granges object do not match."))
  
  # Intersect genomic region and surveyed region
  if (is.null(cnv)) {
    surveyed_region = intersect(surveyed, region, ignore.strand = TRUE)
    surveyed_region_length = sum(width(surveyed_region))
  } else {
    intersect = intersect(surveyed, region, ignore.strand = TRUE)
    ov = findOverlaps(surveyed, intersect, ignore.strand = TRUE) # query, subject
    surv_overlap = surveyed[queryHits(ov)] 
    # queryHits and subjectHits always return same number of indexes
    # smart in case of multiple hits
    surveyed_region = intersect[subjectHits(ov)]
    mcols(surveyed_region) <- mcols(surv_overlap)
    #https://bioinformatics.stackexchange.com/questions/874/intersection-of-two-genomic-ranges-to-keep-metadata
    surveyed_region_length = sum(as.numeric(width(surveyed_region)) * as.numeric(mcols(surveyed_region)$cn))
  }
  
  # Find which mutations lie in surveyed genomic region
  overlap = findOverlaps(vcf, surveyed_region)
  muts_in_region = as.data.frame(as.matrix(overlap))$queryHits
  
  observed = length(muts_in_region)
  prob = n_muts / surveyed_length
  expected = prob * surveyed_region_length
  
  res = data.frame(n_muts,
                   surveyed_length,
                   prob, surveyed_region_length,
                   expected,
                   observed)
  return(res)
}


if (grepl(',', input_bed)) {
  beds <- strsplit(input_bed, ',')
  regions_list <- lapply(beds, import)
  names(regions_list) <- unlist(strsplit(input_bed_name,','))
} else {
  regions <- import(input_bed)
  regions_list <- GRangesList(regions)
  names(regions_list) <- c(input_bed_name)
}

surveyed <- import(input_covered)
surveyed_list <- rep(list(surveyed), length(vcfs))
# warnings on different chr we don't care

# todo load list in the order dictated by vcfs
cnvs <- list()
wanted_cnvs <- c("1","2","3")
for (name in whole_names) {
  c <- read.table(paste0(input_cnv_dir, "/", name ,".callable.bed.gz"), sep="\t", header=FALSE)
  f <- tempfile()
  write.table(c[,c(1,2,3,4,5)], f, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  c <- import.bed(f) # veeerybad
  # no get the matched ones and remove those outside of cnv 1,2,3
  names(mcols(c)) <- c("cn","cn_other")
  c <- c[mcols(c)$cn %in% wanted_cnvs]
  cnvs[[name]] <- c
}

distr <- my_genomic_distribution(vcfs, surveyed_list, regions_list)
cnv_distr <- my_genomic_distribution(vcfs, surveyed_list, regions_list, cnvs)

#distr_testwhole <- enrichment_depletion_test(distr, by = rep("whole", length(vcfs)))
distr_test <- enrichment_depletion_test(distr, by=sample_names)

pdf(output_plot)
p1 <- plot_enrichment_depletion(distr_test)
#plot_enrichment_depletion(distr_testwhole)

#c_distr_testwhole <- enrichment_depletion_test(cnv_distr, by = rep("whole", length(vcfs)))
c_distr_test <- enrichment_depletion_test(cnv_distr, by=sample_names)
p2 <- plot_enrichment_depletion(c_distr_test)

distr_test$type <- "no_corr_cn"
c_distr_test$type <- "corr_cn"

res <- rbind(distr_test, c_distr_test)
write.table(res, file=output_tsv, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
save.image(paste0(output_tsv,".Rdata"))

grid.arrange(p1, p2)
ggsave(output_plot)
#plot_enrichment_depletion(c_distr_testwhole)

