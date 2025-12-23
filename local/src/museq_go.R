library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(GO.db)
setwd('/scratch/trcanmed/AF_spectra/local/share/data/MUseq')
f <- 'CRC1307.hg38_multianno.txt'

d <- read.table(f, sep='\t', header=T, stringsAsFactors = F)

# get universe
uni <- read.table('/mnt/trcanmed/snaketree/task/variant_annotations/dataset/annovar/hg38/humandb/hg38_refGene.txt', sep='\t', header=F, stringsAsFactors = F)
geneUni <- unique(uni$V13)
# get exonic non syn
#uni <- unique(unlist(strsplit(d$Gene.refGene, ';')))

nonsyn <- d[d$Func.refGene=="exonic" & d$ExonicFunc.refGene=="nonsynonymous SNV",]
 
geneList <- unique(unlist(strsplit(nonsyn$Gene.refGene, ';')))

### as.character universe
geneUni <- gene_univ_df$V1
geneUni <- as.character(geneUni)

#m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
#dplyr::select(gs_name, human_gene_symbol)

#egomf <- enrichGO(gene          = geneList,
#                  universe      = geneUni,
#                  OrgDb         = "org.Hs.eg.db",
#                  keyType = "SYMBOL",
#                  ont           = "BP",
#                  pAdjustMethod = "BH",  
#                  pvalueCutoff  = 1,
#                  qvalueCutoff  = 1,
#                  readable      = FALSE)
#egomf[egomf$p.adjust < 0.05,]

# cell cycle GO:0007049
# apoptosis (cell death) GO:0008219
# dna repair GO:0006281
w <- c('GO:0007049', 'GO:0008219', 'GO:0006281') # all bp

###

go_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = geneList,
  keytype  = "SYMBOL",
  columns  = c("GO", "ONTOLOGY")
)

go_annot$TERM <- mapIds(
  GO.db,
  keys     = go_annot$GO,
  column   = "TERM",
  keytype  = "GOID",
  multiVals = function(x) {paste0(x, collapse=',')}
)

go_annot[go_annot$GO %in% w, ]
#######################################

go_collapsed <- go_annot %>%
  group_by(SYMBOL, ONTOLOGY) %>%
  summarise(
    GO_terms = paste(unique(TERM), collapse = "; "),
    .groups = "drop"
  )

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


results <- genes_hit_go_terms(geneList, w)
