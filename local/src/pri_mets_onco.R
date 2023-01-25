muts <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/mutect_paired/merged.table_nomultiallele_alltiers', header=TRUE, sep="\t", row.names=1)
annot <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/mutect_paired/merged.table_nomultiallele_alltiers_annot',  header=TRUE, sep="\t", row.names=1, stringsAsFactors = F)

table(unique(sapply(annot$cds, function(x){strsplit(x, ":")[[1]][1]})))

all_genes <- unique(annot$genes)

get_genes_binary <- function(sample, AF, genes, thr, all_genes) {
  keep <- AF[,sample] > thr
  genes <- genes[keep,]
  return(all_genes %in% genes$genes)
}
thr <- 0.2
all <- sapply(colnames(muts), get_genes_binary, muts, annot, thr, all_genes)
rownames(all) <- all_genes
pri <- all[,grepl('PRX', colnames(all))]
met <- all[,grepl('LMX', colnames(all))]

tpri <- t(pri)
tmet <- t(met)
colnames(tmet) <- substr(colnames(tmet), 0, 7)
colnames(tpri) <- substr(colnames(tpri), 0, 7)

both <- tmet & tpri
both2 <- t(apply(both, 1, as.numeric))
rownames(both2) <- rownames(both)
p <- ifelse(tpri-both ==1, T, F)
m <- ifelse(tmet-both==1,T, F)

mat_list <- list(Both=t(both), PRX=t(p) , LMX=t(m))


col = c("Both" = "darkgoldenrod3", "PRX" = "darkblue", "LMX"= "firebrick1")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Both = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Both"], col = NA))
  },
  PRX = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["PRX"], col = NA))
  },
  LMX = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["LMX"], col = NA))
  }
)

column_title = "PriMet OncoPrint"


rs <- mat_list[[1]]+mat_list[[2]]+mat_list[[3]]
rsu <- rowSums(rs)
rsu <- rsu[order(-rsu)]

keep <- names(head(rsu, n=50))

mat_list_2 <- lapply(mat_list, function(x, keep) { x[rownames(x) %in% keep,] }, keep)

mat_list_2 <- lapply(mat_list_2, function(x) { colnames(x) <- substr(colnames(x),0,7); return(x) })

s <- mat_list_2[[1]]+mat_list_2[[2]]+mat_list_2[[3]]
su <- colSums(s)
su <- su[order(-su)]

op <- oncoPrint(mat_list_2, alter_fun = alter_fun, col = col, column_order = names(su),
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, pct_gp=gpar(fontsize=5), column_names_gp=gpar(fontsize=5),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=5)))

op

s <- mat_list[[1]]+mat_list[[2]]+mat_list[[3]]
su <- colSums(s)
su <- su[order(-su)]

op <- oncoPrint(mat_list, alter_fun = alter_fun, col = col, column_order = names(su),
                remove_empty_columns = TRUE, remove_empty_rows = TRUE, pct_gp=gpar(fontsize=5), column_names_gp=gpar(fontsize=5),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=5)))

op

