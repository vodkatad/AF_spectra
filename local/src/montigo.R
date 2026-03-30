library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(GO.db)
library(clusterProfiler)
library(msigdbr)
#load('/scratch/trcanmed/AF_spectra/dataset_museq/subclonal/pipGO.Rdata')

set.seed(42)

n <- as.numeric(snakemake@params[['n']])
rep <- as.numeric(snakemake@params[['rep']])
outf <- snakemake@output[['monti']]
print(n)
print(rep)
print(outf)

load(snakemake@input[['pipGO']])

genes_hit_go_terms <- function(gene_vector, w, orgdb = org.Hs.eg.db) {
  
  # Step 1: Convert gene symbols to Entrez IDs
  gene_df <- bitr(
    gene_vector,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = orgdb
  )
  
  # Step 2: Get direct GO annotations
  go_annot <- AnnotationDbi::select(
    orgdb,
    keys     = gene_df$ENTREZID,
    keytype  = "ENTREZID",
    columns  = c("GOALL", "ONTOLOGYALL")
  ) %>%
    left_join(gene_df, by = "ENTREZID")
  
  go_annot <- go_annot[go_annot$GOALL %in% w,]
  return(length(unique(go_annot$SYMBOL)))
}

w <- c('GO:0007049', 'GO:0008219', 'GO:0006281') # all bp
dfg <- data.frame(goid=w, term=c('cell cycle', 'cell death', 'dna repair'))

rand <- function() {
  geneList <- sample(geneUni, n)
  results <- genes_hit_go_terms(geneList, w)
}

monti <- replicate(rep, rand())

saveRDS(monti, file=outf)
