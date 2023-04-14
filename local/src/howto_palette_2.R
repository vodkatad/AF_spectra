
colors <- '/scratch/trcanmed/AFspectra/local/share/data/palette.rds'
palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model_clone
pp <- as.data.frame(table(pal))
pp <- pp[pp$Freq>1,]
palette[palette %in% as.character(pp$pal)]
pal[names(pal)=='CRC1599LM_01']='#ffad33'
pal[names(pal)=='CRC1078_07']="#2fc603"
palette_df <- data.frame(palette=as.character(pal), model_clone=names(pal), stringsAsFactors = FALSE)

write_rds(palette_df, path=colors)
