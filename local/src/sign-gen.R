library(ggplot2)

load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_10.Rdata')

palette_df <- readRDS('/scratch/trcanmed/AF_spectra/local/share/data/model_palette.rds')
pal <- palette_df$palette
names(pal) <- palette_df$model

sign_f <- '/scratch/trcanmed/AF_spectra/datasetV2/vitrovivobulk_heatmap_merged_cosmic.tsv'
gen_f <- '/scratch/trcanmed/AF_spectra/local/share/data/SourceData/gens.tsv'


#todo also PR!
sign <- read.table(sign_f, sep='\t', quote='', stringsAsFactors = F, header=T)
gen <- read.table(gen_f, sep='\t', quote='', stringsAsFactors = F, header=F)

sign <- sign[grepl('vitroMA', rownames(sign)) & !grepl('2nd', rownames(sign)),]
rownames(sign) <- gsub('_vitroMA', '', rownames(sign))

colnames(gen) <- c('model', 'clone', 'gen', 'tipo', 'days')
# remove 2nd round for now
gen <- gen[!grepl('[A-Z]', gen$clone),]
gen <- gen[gen$tipo== 'edu',]
#they are all equal does not matter
avedays <- tapply(gen$days, gen$model, mean)
avegen <- tapply(gen$gen, gen$model, mean)
all(names(avedays)==names(avegen))  
dg <- data.frame(rownames=names(avegen), days=avedays, gen=avegen)

ggplot(data=dg, aes(x=days, y=gen, color=model))+geom_point()+theme_bw(base_size = 20)+
scale_color_manual(values=pal)+ theme_bw(base_size = 12)

colnames(sign) <- paste0('SBS', seq(1,ncol(sign)))
d <- merge(dg, sign, by='row.names')

ggplot(data=d, aes(x=days, y=SBS8,color=model))+geom_point()+
scale_color_manual(values=pal)+ theme_bw(base_size = 12)+theme(legend.position="none")

ggplot(data=d, aes(x=gen, y=SBS8,color=model))+geom_point()+
  scale_color_manual(values=pal)+ theme_bw(base_size = 12)+theme(legend.position="none")


ggplot(data=d, aes(x=days, y=SBS18,color=model))+geom_point()+
  scale_color_manual(values=pal)+ theme_bw(base_size = 12)+theme(legend.position="none")

ggplot(data=d, aes(x=gen, y=SBS18,color=model))+geom_point()+
  scale_color_manual(values=pal)+ theme_bw(base_size = 12)+theme(legend.position="none")


cor.test(d$gen, d$SBS8, method='spearman')
cor.test(d$days, d$SBS8, method='spearman')


cor.test(d$gen, d$SBS18, method='spearman')
cor.test(d$days, d$SBS18, method='spearman')
