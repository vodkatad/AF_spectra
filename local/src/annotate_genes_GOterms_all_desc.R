library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(GO.db)
library(clusterProfiler)


annovar_f <- snakemake@input[["annovar"]]
out_f <- snakemake@output[["outf"]]
log_f <- snakemake@log[['log']]

d <- read.table(annovar_f, sep='\t', header=T, stringsAsFactors = F)
nonsyn <- d[d$Func.refGene=="exonic" & d$ExonicFunc.refGene=="nonsynonymous SNV",]
geneList <- unique(unlist(strsplit(nonsyn$Gene.refGene, ';')))

sink(log_f)
nrow(nonsyn)
geneList
sink()
# cell cycle GO:0007049
# apoptosis (cell death) GO:0008219
# dna repair GO:0006281
w <- c('GO:0007049', 'GO:0008219', 'GO:0006281') # all bp

### downward traversing
# Function: check if genes are annotated to GO terms or any of their descendants
genes_hit_go_terms <- function(gene_vector, go_terms_of_interest, ontology="BP", orgdb = org.Hs.eg.db) {
  
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
    columns  = c("GO", "ONTOLOGY")
  ) %>%
    left_join(gene_df, by = "ENTREZID")
  
  # Step 3: Helper to get direct children of a GO term
  get_go_children <- function(go_id, ontology) {
    if (ontology == "BP") {
      as.list(GOBPCHILDREN)[[go_id]]
    } else if (ontology == "MF") {
      as.list(GOMFCHILDREN)[[go_id]]
    } else if (ontology == "CC") {
      as.list(GOCCCHILDREN)[[go_id]]
    } else character(0)
  }
  
  # Step 4: Recursive function to get all descendants
  # BROKEN FIXME
  get_all_descendants <- function(go_id, ontology, seen = character()) {
    #if (go_id %in% seen) return(character())
    children <- get_go_children(go_id, ontology)
    if (is.na(children[1])) {
      return(seen)
    } else {
      descendants <- unlist(
        lapply(children, get_all_descendants, ontology = ontology, seen = children)
      )
      return(unique(c(children, descendants)))
    }
  }
  
  get_all_descendantsAI <- function(go_id, ontology, seen = character()) {
    if (go_id %in% seen) return(character())
    children <- get_go_children(go_id, ontology)
    descendants <- unlist(
      lapply(children, get_all_descendantsAI, ontology = ontology, seen = children)
    )
    unique(c(children, descendants))
  }
  
  
  # Step 5: Build descendant map for all GO terms of interest
  
  
  go_desc_map <- lapply(go_terms_of_interest, function(go) {
    unique(c(go, get_all_descendantsAI(go_id=go, ontology=ontology)))
  })
  names(go_desc_map) <- go_terms_of_interest
  
  # Step 6: Check which genes hit which terms (or descendants)
  gene_hits <- go_annot %>%
    rowwise() %>%
    mutate(
      HIT_TERM = paste(
        names(go_desc_map)[sapply(go_desc_map, function(desc) GO %in% desc)],
        collapse = ";"
      )
    ) %>%
    ungroup() %>%
    filter(HIT_TERM != "")
  
  return(gene_hits)
}

# =========================
# Example usage:


results <- as.data.frame(genes_hit_go_terms(geneList, w))

write.table(results, file=out_f, sep="\t", row.names=F, quote=F)
