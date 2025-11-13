all <- read.table('/scratch/trcanmed/AF_spectra/dataset_MAtreats2/all_muts.tsv', sep="\t", header=F, stringsAsFactors = F)
colnames(all) <- c('n', 'model')
coding <- read.table('/scratch/trcanmed/AF_spectra/dataset_MAtreats2/all_muts_coding.tsv', sep="\t", header=T, stringsAsFactors = F)

ncoding <- sapply(unique(coding$model), function(x) { m <- coding[coding$model==x,'Freq']; sum(m)})
ncod <- data.frame(model=unique(coding$model), n=ncoding)

m <- merge(ncod, all, by='model')
colnames(m) <- c('model', 'exonic', 'tot')
m$not_exonic <- m$tot - m$exonic

library(reshape)
longf <- melt(m, id='model')
longf <- longf[longf$variable != "tot",]

library(dplyr)
longfp <- longf %>% 
  group_by(model) %>% 
  mutate(perc = value/sum(value))

longfp$tmodel <- as.character(longfp$model)
longfp$model <- sapply(longfp$tmodel, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})
longfp$treat <- sapply(longfp$tmodel, function(x) {y<-strsplit(x, '_')[[1]][1]; return(y[1])})

longfp <- longfp[order(longfp$model, longfp$treat, longfp$treat),]
longfp$order <- seq(1, nrow(longfp))

load('/scratch/trcanmed/AF_spectra/dataset_MAtreats/themeok_6.Rdata')
y_breaks <- guess_ticks(c(0,1), nticks=5)

p <- ggplot(as.data.frame(longfp), aes(x = reorder(tmodel, order), y = perc, fill = variable)) +
  geom_col(color = "black") +
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("mediumspringgreen", "lightgoldenrod4"))+
  unmute_theme + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+
  scale_y_continuous(breaks=y_breaks, limits=c(0,1), expand = c(0, 0))+
  scale_x_discrete(labels=c('NT_CRC1430'="NT", 'T_CRC1430'="Afatinib", 'NT_CRC1620'="NT", 'T_CRC1620'="Afatinib"))+
  xlab('')+ylab('% accumulated SNVs')


coding <- coding[!grepl('baseline', coding$model),]
longfp <- coding %>% 
  group_by(model) %>% 
  mutate(perc = Freq/sum(Freq))

longfp$tmodel <- as.character(longfp$model)
longfp$model <- sapply(longfp$tmodel, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})
longfp$treat <- sapply(longfp$tmodel, function(x) {y<-strsplit(x, '_')[[1]][1]; return(y[1])})

longfp <- longfp[order(longfp$model, longfp$treat, longfp$treat),]
longfp$order <- seq(1, nrow(longfp))

y_breaks <- guess_ticks(c(0,1), nticks=5)

p <- ggplot(as.data.frame(longfp), aes(x = reorder(tmodel, order), y = perc, fill = Var1)) +
  geom_col(color = "black") +
  geom_text(aes(label = Freq),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('lightgray', "aquamarine", "steelblue", 'slategray4'))+
  unmute_theme + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))+
  scale_y_continuous(breaks=y_breaks, limits=c(0,1), expand = c(0, 0))+
  scale_x_discrete(labels=c('NT_CRC1430'="NT", 'T_CRC1430'="Afatinib", 'NT_CRC1620'="NT", 'T_CRC1620'="Afatinib"))+
  xlab('')+ylab('% accumulated SNVs')


