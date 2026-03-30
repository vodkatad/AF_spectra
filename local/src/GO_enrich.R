library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(GO.db)
library(clusterProfiler)
library(msigdbr)

genes_f <- snakemake@input[["genes"]]
universe_f <- snakemake@input[["universe"]]
enrich_f <- snakemake@output[["enrich"]]

d <- read.table(genes_f, sep='\t', header=T, stringsAsFactors = F)
geneList <- unique(as.character(d$V1))

## standard enrichments
# get exonic non syn
#uni <- unique(unlist(strsplit(d$Gene.refGene, ';')))

### as.character universe
uni <- read.table(universe_f, sep='\t', header=F, stringsAsFactors = F)
geneUni <- unique(uni$V1)
geneUni <- as.character(geneUni)

#m_t2g <- msigdbr(collection = "C6", db_species='HS', species='human') %>% 
#dplyr::select(gs_name, gene_symbol)

ego <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "ALL",
                pAdjustMethod = "bonferroni",  
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)

if (!is.null(ego)) {
  # readjust cause it does each ontology separately with ALL
  ego <- as.data.frame(ego)
  ego$padj <- p.adjust(ego$pvalue, method='bonferroni')
  res <- ego[ego$padj < 0.05,]
  write.table(res, file=enrich_f, sep='\t', quote=F, row.names=FALSE)
} else {
  sink(log_f)
  'No GO'
  sink()
  write.table(data.frame(), file=enrich_f, sep='\t', quote=F, row.names=FALSE)
}

##
