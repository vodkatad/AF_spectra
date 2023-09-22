library(ggplot2)

#setwd('/home/data/Dropbox/work/evol/MA/paired_patients_nQ/fits')
dnobin <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/bestbet_0.12_0.24.tsv', sep="\t", header=TRUE)
#d <- read.table('bestbet_bin_0.12_0.24.tsv', sep="\t", header=TRUE)

#d$lmodel <- substr(rownames(d), 0, 10)
dnobin$lmodel <- substr(rownames(dnobin), 0, 10)

#mm <- merge(d, dnobin, by="lmodel")

d <- dnobin

d$smodel <- substr(rownames(d), 0, 7)
d$mp <- substr(rownames(d), 8, 10)

d$mp <- factor(d$mp, levels=c('PRX', 'LMX'))

d2 <- d[d$r > 0.90 & d$subcl > 10,]

dd <- as.data.frame(table(d2$smodel))
d2 <- d2[d2$smodel %in% dd[dd$Freq == 2,'Var1'],]

ggplot(data=d2, aes(x=smodel, y=intercept, fill=mp))+
  geom_col(position="dodge")+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Patient')+ylab('Slope μ/β')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))

lar <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'intercept'] - sub[sub$mp=="PRX",'intercept']
}


meslo <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'intercept']
}


dd2 <- as.data.frame(sapply(unique(d2$smodel), lar, d2))
dd2 <- dd2[order(dd2[,1]),, drop=FALSE]

dd3 <- as.data.frame(sapply(unique(d2$smodel), meslo, d2))
dd3 <- dd3[order(dd3[,1]),, drop=FALSE]



d3 <- merge(dd2, d2, by.x="row.names", by.y="smodel")
colnames(d3)[2] <- 'o'
colnames(d3)[1] <- 'smodel'
ggplot(data=d3, aes(x=reorder(smodel,o), y=intercept, fill=mp))+
  geom_col(position="dodge")+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Patient')+ylab('Slope μ/β')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))

d <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/univT0_pileup/all.MR_ov', sep="\t", header=T)
d2 <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/CRC1307/univMutect/all.MR_ov', sep="\t", header=T)
m <- merge(d, d2, by.x="start", by.y="end")
m$MR_EDUpile <- m$MR_EDU.x * 1000000000
m$MR_EDUmute <- m$MR_EDU.y * 1000000000

ggplot(data=m,aes(x=MR_EDUmute, y=MR_EDUpile))+geom_point()+theme_bw()+xlab('univMutect')+ylab('univpileup')+geom_abline(slope=1, intercept=0)+xlim(1.25,1.75)+ylim(1.25,1.75)

#### New figure 5a hypothetical
## cumulative fits
load('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/CRC1599PRX0A02002TUMD03000V2.fit.0.12_0.24.pdf.debug.RData')
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", ylim=c(0,65))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model, col="#ffcc33")
print(sfit$r.squared)

invf_p <- invf
excum_p <- excum
model_p <- model
load('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/CRC1599LMX0A02001TUMD03000V2.fit.0.12_0.24.pdf.debug.RData')
plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", ylim=c(0,65))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model, col="#ff9900")
print(sfit$r.squared)


data <- data.frame(n=c(excum_p, excum), invf=c(invf_p, invf), type=c(rep('PRX', length(excum_p)), rep('LMX', length(excum))))
x_breaks<-guess_ticks(data$invf)
y_breaks<-guess_ticks(data$n)

data$type <- factor(data$type, levels= c('PRX', 'LMX'))

ggplot(data=data, aes(x=invf, y=n))+geom_point()+theme_bw()+stat_smooth(aes(color=type), method = "lm")+ 
  xlab('1/f')+ ylab("Cumulative n. of muts M(f)")+scale_color_manual(values=c('#adacac', '#595959'))+
  unmute_theme

p <- ggplot(data=data, aes(x=invf, y=n, color=type))+geom_point()+stat_smooth(method = "lm")+ 
  xlab('1/f')+ ylab("Cumulative n. of muts M(f)")+scale_color_manual(values=c('#adacac', '#595959'))+
  theme_bw()+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)), expand = c(0, 0))

ggsave('~/5b_alternate.svg', plot=p, width=60, height=60, units="mm")


load ('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_5c_slopes_pairedscatter.svg.Rdata')
library(reshape)


me <- cast(data=fit_r2, 'smodel~mp', value="intercept")
x_breaks<-guess_ticks(me$PRX, fixed_max=68)
y_breaks<-guess_ticks(me$LMX, fixed_max=68)


p <- ggplot(data=me, aes(x=PRX, y=LMX)) + geom_point(shape=1, size=1) + geom_abline(slope=1, intercept=0)+
scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)), expand = c(0, 0))+
theme(text=element_text(size=15))+unmute_theme

ggsave('~/5c_alternate.svg', plot=p, width=60, height=60, units="mm")
