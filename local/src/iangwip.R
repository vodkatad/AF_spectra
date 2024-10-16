library(ggplot2)
age <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/ages.txt', sep="\t", header=F)
colnames(age) <- c('model', 'age')

sign <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/vitrovivobulk_heatmap_merged_cosmic.tsv', sep="\t", header=T)

palette_df <- readRDS('/scratch/trcanmed/AF_spectra/local/share/data/IANG_palette.rds')
pal <- palette_df$palette
names(pal) <- palette_df$model


sign <- sign[grepl('bulk', rownames(sign)), ]
sign$model <- gsub("_bulk", "", rownames(sign))

sign$sample <- sign$model
sign$model <- paste0(sign$model, ifelse(!grepl('\\d$', sign$model), '', ifelse(sign$model=="CRC0282", 'PR', 'LM')))


#treat <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/ma_treats_cleverers_bulk/all_sign_treat.tsv', header=F, sep="\t")
#colnames(treat) <- c('model', paste0('X', seq(1, 30)))

sign <- sign[, c('model',  'X1')]
treat <- treat[, c('model', 'X1')]

#sign <- rbind(sign, treat)

m <- merge(age, sign, by="model")

#ggplot(data=m ,aes(x=X1, y=age, color=sample))+geom_point()+theme_bw(base_size=15)+xlab('SBS1')+scale_color_manual(values=pal)+theme(legend.position = 'none')


iang <- read.table('/scratch/trcanmed/AF_spectra/dataset_IANG/vitrobulk_heatmap_merged_cosmic.tsv', sep="\t", header=T)


colors <- "#10d8eb,#5b3d80,#e69c6e"
models <- "CRCUECHPR,CRC2826PR,CRC3023PR"
model_palette_df <- data.frame(palette=unlist(strsplit(colors, ',')), model=unlist(strsplit(models, ',')), stringsAsFactors=F)
pal <- model_palette_df$palette
names(pal) <- model_palette_df$model

sub <- iang[grepl('bulk', rownames(iang)),]
sub$name <- gsub('_bulk', '', rownames(sub))
sub$name <- gsub('PRO', 'PR', sub$name)

m$age <- as.numeric(m$age)

fit <- lm(data=m, formula=as.formula(age~X1))
sfit <- summary(fit)
p1 <- ggplot()+geom_point(data=m ,aes(x=X1, y=age))+
  geom_vline(data=sub, aes(xintercept=X1, color=name), lwd=1)+
  theme_bw(base_size=15)+xlab('SBS1')+theme(legend.position='none')+scale_color_manual(values=pal)+geom_smooth(data=m ,aes(x=X1, y=age),method=lm, fullrange=T, lwd=0.5)+
  geom_abline(slope=sfit$coefficients[2,1], intercept=sfit$coefficients[1,1], color='blue', lwd=0.5)


sub$age <- c(23, 50, 50)
p2 <- ggplot()+geom_point(data=m ,aes(x=X1, y=age))+
  geom_point(data=sub, aes(x=X1, y=age, color=name, size=1))+
  theme_bw(base_size=15)+xlab('SBS1')+scale_color_manual(values=pal)+geom_smooth(data=m ,aes(x=X1, y=age),method=lm, fullrange=T, lwd=0.5)+
  geom_abline(slope=sfit$coefficients[2,1], intercept=sfit$coefficients[1,1], color='blue', lwd=0.5)

ggsave(plot=p1, file='/scratch/trcanmed/AF_spectra/dataset_fIANG/lie1.svg',width=89, height=89, units="mm")
ggsave(plot=p2, file='/scratch/trcanmed/AF_spectra/dataset_fIANG/lie2_legend.svg',width=89, height=89, units="mm")


p3 <- ggplot()+geom_point(data=m ,aes(x=X1, y=age))+
  geom_point(data=sub, aes(x=X1, y=age, color=name, size=1))+
  theme_bw(base_size=15)+theme(legend.position='none')+xlab('SBS1')+scale_color_manual(values=pal)+geom_smooth(data=m ,aes(x=X1, y=age),method=lm, fullrange=T, lwd=0.5)+
  geom_abline(slope=sfit$coefficients[2,1], intercept=sfit$coefficients[1,1], color='blue', lwd=0.5)

ggsave(plot=p3, file='/scratch/trcanmed/AF_spectra/dataset_fIANG/lie2.svg',width=89, height=89, units="mm")


colors <- "#10d8eb,#5b3d80,#e69c6e"
models <- "CRCUECHPRO,CRC2826PRO,CRC3023PRO"
model_palette_df <- data.frame(palette=unlist(strsplit(colors, ',')), model=unlist(strsplit(models, ',')), stringsAsFactors=F)
saveRDS(model_palette_df, file="/scratch/trcanmed/AF_spectra/local/share/data/IANG_palette.rds")


models <- c("CRCUECHPRO-01","CRCUECHPRO-05","CRCUECHPRO-06","CRCUECHPRO-07","CRCUECHPRO-13", "CRC2826PRO-03","CRC2826PRO-05","CRC2826PRO-09","CRC2826PRO-14", "CRC3023PRO-02")
colors <- "#10d8eb,#0fc0d1,#097c87,#47dfed,#7ec5cc,#5b3d80,#a772e8,#7e57ad,#643f91,#e69c6e"
model_palette_df <- data.frame(palette=unlist(strsplit(colors, ',')), model_clone=unlist(strsplit(models, ',')), stringsAsFactors=F)
saveRDS(model_palette_df, file="/scratch/trcanmed/AF_spectra/local/share/data/IANG_allclones_palette.rds")

##
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_1b_MR.svg.Rdata')
y_breaks <- guess_ticks(c(0,4))

pl <- ggplot() + 
  geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.7))+
  geom_segment(data=pdata, aes(x=xmodel-0.2, yend=mean,y=mean,  xend=xmodel+0.2),size=.3) +
  geom_errorbar(data=pdata, aes(ymin=lower, ymax=upper, x=model), size=0.3, width=0.3)+ylab('MR [SNV/(Gbp*division)]')+xlab('PDTs')+
  scale_color_manual(values=pal)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   


ggsave(pl, file='~/trimmed.svg')

y_breaks <- guess_ticks(c(0,2))

pl <- ggplot() + 
  geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.7))+
  geom_segment(data=pdata, aes(x=xmodel-0.2, yend=mean,y=mean,  xend=xmodel+0.2),size=.3) +
  geom_errorbar(data=pdata, aes(ymin=lower, ymax=upper, x=model), size=0.3, width=0.3)+ylab('MR [SNV/(Gbp*division)]')+xlab('PDTs')+
  scale_color_manual(values=pal)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   


ggsave(pl, file='~/trimmed2.svg')


ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('Gained muts')+ylab(ylab)+xlab('')+
  geom_point(data=our, aes(x=model, y=gained, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=pal)+scale_shape_manual(values=scale_shape)


ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('Gained muts')+ylab(ylab)+xlab('')+
  geom_point(data=our, aes(x=model, y=gained, color=model_clone), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=pal)

