df <- data.frame(palette=c('#c05666', '#56adc0', "#aec056"), model=c('CRC2608PR', 'CRC2566LM', 'CRC2573LM'), stringsAsFactors=FALSE)
saveRDS(df, file='../local/share/data/cleverers_palette.rds')
savehistory('../local/share/data/cleverers_palette.R')
