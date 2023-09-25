lm_fit  <- snakemake@input[['LM']]
pr_fit  <- snakemake@input[['PR']]

outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))
library(ggplot2)
load(theme)


load(pr_fit)
#plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", ylim=c(0,65))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
print(oi[labels])
print(oex[labels])
#axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
#abline(model, col="#ffcc33")
#print(sfit$r.squared)

#> oex[labels]
#[1] 0.238 0.222 0.192    NA 0.136
#> oi[labels]
#[1] 0.03501401 0.33783784 1.04166667         NA 3.18627451
# TODO INDAGARE perchè? Labels ggplot x axis needs fixing

invf_p <- invf
excum_p <- excum
model_p <- model

load(lm_fit)
#plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", ylim=c(0,65))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
#axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
#abline(model, col="#ff9900")
#print(sfit$r.squared)

data <- data.frame(n=c(excum_p, excum), invf=c(invf_p, invf), type=c(rep('PRX', length(excum_p)), rep('LMX', length(excum))))
x_breaks<-guess_ticks(data$invf)
y_breaks<-guess_ticks(data$n, fixed_max=64)

data$type <- factor(data$type, levels= c('PRX', 'LMX'))

p <- ggplot(data=data, aes(x=invf, y=n, color=type))+geom_point(shape=1, size=1)+stat_smooth(method = "lm", se=FALSE, show.legend = FALSE, size=0.5)+ 
  xlab('1/f')+ ylab("Cumulative n. of muts M(f)")+scale_color_manual(values=c('#ffcc33', '#ff9900'))+
  unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)), expand = c(0, 0))
  
pdf('fig_5b_cumfit.pdf')
print(p)
graphics.off()
ggsave(outplot, plot=p, width=60, height=60, units="mm")


save.image(paste0(outplot, '.Rdata'))