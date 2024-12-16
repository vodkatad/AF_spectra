sign_f <- '/scratch/trcanmed/AF_spectra/datasetV2/vitrovivobulk_heatmap_merged_cosmic.tsv'
MR_f <- '/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/sorted_MR_avg_new.tsv'
mut_bur_f <- '/scratch/trcanmed/AF_spectra/datasetV2/all_clonal_0_n.txt'
  #/local/share/data/bulkburdens_indel

sign <- read.table(sign_f, sep="\t", header=T, stringsAsFactors = F)
MR <- read.table(MR_f, sep="\t", header=T, stringsAsFactors = F)
mut_bur <- read.table(mut_bur_f, sep="\t", header=F, stringsAsFactors = F)

sign <- sign[grepl('bulk', rownames(sign)),]
sign_df <- data.frame(sample=rownames(sign), SBS1=sign$X1, stringsAsFactors = F)
sign_df$model <- sapply(sign_df$sample, function(x) {y<-strsplit(x, '_')[[1]][1]; return(y[1])})

mut_bur <- mut_bur[grepl('bulk', mut_bur$V2),]
mut_bur_df <- data.frame(sample=mut_bur$V2, mut_burden=mut_bur$V1, stringsAsFactors = F)
mut_bur_df$model <- sapply(mut_bur_df$sample, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})

m <- merge(sign_df, MR, by="model")
m2 <- merge(m, mut_bur_df, by="model")
m2$estiMR <- m2$mut_burden / m2$SBS1


ggplot(data=m2, aes(x=estiMR, y=mean))+geom_point()

load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_4c_cor.svg.Rdata')
m2$MR<-m2$mean*1000000000

m2$model <- paste0(m2$model, ifelse(!grepl('\\d$', m2$model), '', ifelse(m2$model=="CRC0282", 'PR', 'LM')))

ci <- cor.test(m2$estiMR, m2$MR, method="spearman")


y_breaks<-guess_ticks(m2$MR)
#x_breaks<-guess_ticks(m2$estiMR)

ggplot(m2, aes(x=estiMR, y=MR)) +  geom_point(aes(color=model), size=2)+
  unmute_theme+labs(caption=paste0('spearman=', round(ci$estimate,2), ' pval=',round(ci$p.value, 4))) + scale_color_manual(values=pal)+xlab('mutbur/SBS1')+ylab('MR')+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+theme_bw(base_size=15)