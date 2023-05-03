df <- data.frame(palette=c('#c00000', '#0011c0', "#00e415"), model=c('P1', 'P2', 'P3'), stringsAsFactors=FALSE)
saveRDS(df, file='../local/share/data/clevers_palette.rds')
savehistory('../local/share/data/clevers_palette.R')