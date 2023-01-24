mb <- read.table('/home/egrassi/MA_burdens_01.tsv', sep="\t", header=F)
colnames(mb) <- c('id', 'tot', 'mb')
mb$smodel <- substr(mb$id, 0, 7)
mb[mb$id=="CRC1599PRO-0-B", 'smodel'] <- 'CRC1599PR'
mb[mb$id=="CRC1599LMO-0-B", 'smodel'] <- 'CRC1599LM'
d <- read.table('/scratch/trcanmed/AF_spectra/dataset/MR_edu_SNV_averaged.tsv', sep="\t", header=T)

COLORS_MODELS_NOMSI="#f607b9,#9900ff,#155d00,#77a003,#0829fc,#ff9900,#ffff00" # without 2nd round, aggregated
COLORS_MODELS_2="#cc3300,#f607b9,#9900ff,#155d00,#77a003,#0829fc,#ff9900,#ffff00" # without 2nd round, aggregated

colors <- COLORS_MODELS_2

cbPalette <- unlist(strsplit(colors, ','))

m <- merge(mb, d, by.x="smodel", by.y="model")

ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15), 
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15),  axis.title.x=element_text(size=20),
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

name_x <- 'mean'
name_y <- 'mb'
ci <- cor.test(m[, name_x], m[, name_y])



ggplot(m, aes_string(x=name_x, y=name_y)) +  geom_point(aes(color=smodel), size=3) + geom_smooth(method='lm')+
  ctheme+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=cbPalette)

m <- m[m$smodel !="CRC0282",]

colors <- COLORS_MODELS_NOMSI

cbPalette <- unlist(strsplit(colors, ','))
ci <- cor.test(m[, name_x], m[, name_y])

ggplot(m, aes_string(x=name_x, y=name_y)) +  geom_point(aes(color=smodel), size=3) + geom_smooth(method='lm')+
  ctheme+labs(caption=paste0('pearson=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=cbPalette)
