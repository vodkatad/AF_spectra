library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(GO.db)
library(clusterProfiler)
library(msigdbr)

annovar_f <- snakemake@input[["annovar"]]
universe_f <- snakemake@input[["universe"]]
out_f <- snakemake@output[["outf"]]
enrich_f <- snakemake@output[["enrich"]]
log_f <- snakemake@log[['log']]

d <- read.table(annovar_f, sep='\t', header=T, stringsAsFactors = F)
nonsyn <- d[d$Func.refGene=="exonic" & d$ExonicFunc.refGene=="nonsynonymous SNV",]
geneList <- unique(unlist(strsplit(nonsyn$Gene.refGene, ';')))

## standard enrichments
uni <- read.table('/mnt/trcanmed/snaketree/task/variant_annotations/dataset/annovar/hg38/humandb/hg38_refGene.txt', sep='\t', header=F, stringsAsFactors = F)
geneUni <- unique(uni$V13)
# get exonic non syn
#uni <- unique(unlist(strsplit(d$Gene.refGene, ';')))

### as.character universe
uni <- read.table(universe_f, sep='\t', header=F, stringsAsFactors = F)
geneUni <- unique(uni$V13)
geneUni <- as.character(geneUni)

#m_t2g <- msigdbr(collection = "C6", db_species='HS', species='human') %>% 
#dplyr::select(gs_name, gene_symbol)

ego <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "ALL",
                pAdjustMethod = "BH",  
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)

if (!is.null(ego)) {
  # readjust cause it does each ontology separately with ALL
  ego <- as.data.frame(ego)
  ego$padj <- p.adjust(ego$pvalue, method='BH')
  res <- ego[ego$padj < 0.05,]
  write.table(res, file=enrich_f, sep='\t', quote=F, row.names=FALSE)
} else {
  sink(log_f)
  'No GO'
  sink()
  write.table(data.frame(), file=enrich_f, sep='\t', quote=F, row.names=FALSE)
}

save.image('pipGO.Rdata')

##

sink(log_f, append=TRUE)
nrow(nonsyn)
geneList
sink()
# cell cycle GO:0007049
# apoptosis (cell death) GO:0008219
# dna repair GO:0006281
w <- c('GO:0007049', 'GO:0008219', 'GO:0006281') # all bp
dfg <- data.frame(goid=w, term=c('cell cycle', 'cell death', 'dna repair'))

### downward traversing
# Function: check if genes are annotated to GO terms or any of their descendants
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
  return(go_annot)
}

# =========================
# Example usage:


w <- c('GO:0007049', 'GO:0008219', 'GO:0006281') # all bp
dfg <- data.frame(goid=w, term=c('cell cycle', 'cell death', 'dna repair'))

results <- as.data.frame(genes_hit_go_terms(geneList, w))
results <- merge(results, dfg, by.x='GOALL', by.y='goid')

write.table(results, file=out_f, sep="\t", row.names=F, quote=F)