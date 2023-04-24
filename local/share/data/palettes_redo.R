old_model <- readRDS('old_pal/model_palette.rds')
model <- readRDS('old_pal/palette.rds')
old_model[old_model$model=="CRC1599PR", 'palette'] = '#ffcc33' # '#e8aa0d'
saveRDS(old_model, 'model_palette.rds')
model[model$model=="CRC1599PR_01", 'palette'] = '#ffcc33'
model[model$model=="CRC1599PR_10", 'palette'] = '#cc9900'
saveRDS(model, 'palette.rds')
savehistory('palettes_redo.R')
