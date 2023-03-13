library(ggplot2)
subcl <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/fixedthr_subclonal', header=FALSE, sep="\t")
mredu <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/MR_edu_SNV', header=FALSE, sep="\t")

colnames(subcl) <- c('clone', 'n_subcl')
colnames(mredu) <- c('clone', 'MR')

colors <- '/scratch/trcanmed/AF_spectra/local/share/data/palette.rds'
palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model_clone


m <- merge(subcl, mredu, by= "clone")
m$clone <- as.character(m$clone)
m$cl <- ifelse(grepl('1599', m$clone), substr(m$clone, 1, 12), substr(m$clone, 1, 10))
m$cl <- gsub('-','_', m$cl)
ggplot(data=m, aes(x=MR, y=n_subcl, color=cl))+geom_point()+scale_color_manual(values=pal)+theme_bw()
