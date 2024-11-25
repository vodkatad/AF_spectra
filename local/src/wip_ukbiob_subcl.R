library(ggplot2)
load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/ukbiobank_subcl_cl.svg.Rdata')

table(d$Signature.8_cl == 0)
table(d$Signature.8_subcl == 0)

d$r_cl <- d$Signature.8_cl/d$Signature.1_cl
d$r_sl <- d$Signature.8_subcl/d$Signature.1_subcl


d_nona <- d[!is.infinite(d$r_cl) & !is.infinite(d$r_sl) & !is.na(d$r_cl) & !is.na(d$r_sl),]
pd <- data.frame(ratio=c(d_nona$r_cl, d_nona$r_sl), class=c(rep('Clonal', nrow(d_nona)), rep('Subclonal', nrow(d_nona))))

y_breaks <- guess_ticks(pd$ratio)

p <- ggplot(data=pd, aes(y=ratio, x=class))+geom_boxplot(outlier.size=1)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  unmute_theme+ylab('SBS8/SBS1')+xlab('SNVs')+theme_bw(base_size=20)
##
w <- wilcox.test(d$ratio_subcl, d$ratio_cl, alternative="greater", paired=TRUE)
tt <- wilcox.test(d$lr_sl, d$lr_cl, alternative="greater", paired=TRUE)

d$lr_cl <- log((d$Signature.8_cl/d$Signature.1_cl)+1)
d$lr_sl <- log((d$Signature.8_subcl/d$Signature.1_subcl)+1) 

d_nona <- d[!is.infinite(d$lr_cl) & !is.infinite(d$lr_sl) & !is.na(d$lr_cl) & !is.na(d$lr_sl),]

pd2 <- data.frame(ratio=c(d_nona$lr_cl, d_nona$lr_sl), class=c(rep('Clonal', nrow(d_nona)), rep('Subclonal', nrow(d_nona))))


y_breaks <- guess_ticks(pd2$ratio)

p2 <- ggplot(data=pd2, aes(y=ratio, x=class))+geom_boxplot(outlier.size=1)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+
  unmute_theme+ylab('log(SBS8/SBS1)+1)')+xlab('SNVs')+theme_bw(base_size=20)

##

d$lr_cl <- log((d$Signature.8_cl/d$Signature.1_cl)+0.001)
d$lr_sl <- log((d$Signature.8_subcl/d$Signature.1_subcl)+0.001) 

d_nona <- d[!is.infinite(d$lr_cl) & !is.infinite(d$lr_sl) & !is.na(d$lr_cl) & !is.na(d$lr_sl),]

pd3 <- data.frame(ratio=c(d_nona$lr_cl, d_nona$lr_sl), class=c(rep('Clonal', nrow(d_nona)), rep('Subclonal', nrow(d_nona))))


y_breaks <- guess_ticks(pd2$ratio)

p3 <- ggplot(data=pd3, aes(y=ratio, x=class))+geom_boxplot(outlier.size=1)+
  unmute_theme+ylab('log(SBS8/SBS1)+0.001)')+xlab('SNVs')+theme_bw(base_size=20)
