p1 <- readRDS('../local/share/data/cleverers_palette.rds')
p2 <- readRDS('../local/share/data/clevers_palette.rds')
p1
p2
p1$palette <- p2$palette
#saveRDS(p1, file='../local/share/data/cleverers_palette.rds')
saveRDS(p1, file='../local/share/data/cleverers_palette.rds')
savehistory('../local/share/data//notpastel_cleverers_palette.R')
