library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(GO.db)
library(clusterProfiler)
library(msigdbr)

out_f <- snakemake@output[["outf"]]

w <- c('GO:0007049', 'GO:0008219', 'GO:0006281') # all bp
dfg <- data.frame(goid=w, term=c('cell cycle', 'cell death', 'dna repair'))

go_terms_of_interest <- w
ontology <- "BP"
orgdb <- org.Hs.eg.db
  

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
get_all_descendants <- function(go_id, ontology, seen = character()) {
  if (go_id %in% seen) return(character())
  
  seen <- c(seen, go_id)
  children <- get_go_children(go_id, ontology)
  
  if (length(children) == 0 || is.null(children) || is.na(children[1])) return(character())
  
  children <- as.character(children)
  
  descendants <- unlist(
    lapply(children, function(child) {
      get_all_descendants(child, ontology = ontology, seen = seen)
    }),
    use.names = FALSE
  )
  
  unique(c(children, descendants))
}

# Step 5: Build descendant map for all GO terms of interest


go_desc_map <- lapply(go_terms_of_interest, function(go) {
  unique(c(go, get_all_descendants(go_id=go, ontology=ontology)))
})
names(go_desc_map) <- go_terms_of_interest

# Step 6: Check which genes hit which terms (or descendants)
hit_tbl <- bind_rows(lapply(names(go_desc_map), function(parent_go) {
  tibble(
    GO = go_desc_map[[parent_go]],
    HIT_TERM = parent_go
  )
})) %>%
  distinct()

hit_tbl <- merge(hit_tbl, dfg, by.x='HIT_TERM', by.y='goid')

saveRDS(hit_tbl, file=out_f)