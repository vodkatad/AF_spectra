library(annotatr)
library(rtracklayer)
library(ggplot2)
#bed_gained <- '/scratch/trcanmed/AF_spectra/datasetV2/gained_first.bed'
#bed_bulk <- '/scratch/trcanmed/AF_spectra/datasetV2/gained_bulk.bed'
bed_gained <- '/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/tree/bulk.bed'
bed_bulk <- '/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/platypus_nobin_00/all.gain.bed'

custom_annot <- 'no'
#bed <- snakemake@input[['gained_bed']]
#custom_annot <- snakemake@params[['custom_annot']]
if (custom_annot == "yes") {
  annot_beds <- snakemake@input[['annot_beds']]  
}
output_plot_n <- snakemake@output[['plot_n']]
output_plot_corr <- snakemake@output[['plot_corr']]

if (custom_annot=="no") {
  annots <- c('hg38_basicgenes', 'hg38_genes_intergenic')
  #'hg38_genes_intronexonboundaries')
}  else {
  beds <- read.table(annot_beds, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  for (i in seq(1, nrow(beds))) {
    name <- strsplit(beds[i,1], "_")[[1]][5]
    read_annotations(con = beds[i, 1], genome = 'hg38', name = name, format = 'bed')
  }
  annots <- annotatr_cache$list_env()
  print(annots)
}
annotations <- build_annotations(genome = 'hg38', annotations = annots)
seqlevels(annotations, pruning.mode="coarse") <- paste0('chr', seq(1,22))

my_annotate_regions <- function(bed, annotations) {
  gained <- read_regions(con = bed, genome = 'hg38', format='bed')
  #"coarse": Remove the elements in x where the seqlevels to drop are in use. We should have only chr1-22 anyway, but all the _alt and
  # _random chr from hg38 are still listed as seqlevels and for tidyness we remove them (since we won't evaluate gained muts there does not
  # make sense to possibly randomize them there).
  seqlevels(gained, pruning.mode="coarse") <- paste0('chr', seq(1,22))

  # Intersect the regions we read in with the annotations
  gained_annotated <- annotate_regions(
    regions = gained,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = TRUE)
  
  real_annsum <- summarize_annotations(annotated_regions=gained_annotated, quiet=T)
  
  return(real_annsum)
}

real_annsum <- my_annotate_regions(bed_gained, annotations)
bulk_annsum <- my_annotate_regions(bed_bulk, annotations)

ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.text = element_text(size = 15),
                           legend.title = element_text(size = 15))


real <- as.data.frame(real_annsum)
bulk <- as.data.frame(bulk_annsum)
real$SNV <- 'gained'
bulk$SNV <- 'parental'

pd <- rbind(real, bulk)
pd$annot.type <- as.factor(pd$annot.type)

# I would like to order the annotations from the one with more muts to the least ones, this will be the same in all?
#> levels(pd$annot.type)
#[1] "hg38_genes_1to5kb"               "hg38_genes_3UTRs"                "hg38_genes_5UTRs"               
#[4] "hg38_genes_exons"                "hg38_genes_intergenic"           "hg38_genes_intronexonboundaries"
#[7] "hg38_genes_introns"
names_regions <- levels(pd$annot.type)          
if (custom_annot=="no") {
  renamed_lev <- gsub('hg38_genes_', '', names_regions) # c('1to5kb','3UTRs', '5UTRs','exons','intergenic','introns','promoters')
} else {
  renamed_lev <- gsub('hg38_custom_', '', names_regions) 
}
levels(pd$annot.type) <- renamed_lev
pd$annot.type <- as.character(pd$annot.type) # to unique inheriting the order from n
orderpd <- pd[pd$SNV=="gained",]
orderpd <- orderpd[order(-orderpd$n),]
pd$annot.type <- factor(pd$annot.type, levels=unique(orderpd$annot.type))

ggplot(data=pd, aes(x=annot.type,y=n,fill=SNV))+geom_bar(stat="identity", position='dodge')+ctheme+xlab('')+ylab('#')+
  scale_fill_manual(values=c('darkgreen', 'darkgoldenrod'))


ns <- data.frame(row.names=pd[pd$SNV == "gained", "annot.type"], gained=pd[pd$SNV == "gained", "n"],
                 bulk=pd[pd$SNV == "parental", "n"])

chisq.test(ns)

library(dplyr)

df <- pd %>% 
  group_by(SNV) %>% # Variable to be transformed
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))


ggplot(df, aes(x = "", y = perc, fill = annot.type)) +
  geom_col() +
  coord_polar(theta = "y")+facet_grid(~SNV)+theme_light()


ggplot(df, aes(x = "", y = perc, fill = annot.type)) +
  geom_col() +facet_grid(~SNV)+theme_light()
