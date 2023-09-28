wanted <- c('1', '8', '18')
data <- as.data.frame(data)

wdata <- data[rownames(data) %in% paste0('Signature.', wanted),]

twdata <- t(wdata)
normalized <- twdata/rowSums(twdata)