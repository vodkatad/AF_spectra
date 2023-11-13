lm_fit  <- snakemake@input[['LM']]
pr_fit  <- snakemake@input[['PR']]

outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))
library(ggplot2)
load(theme)


load(pr_fit)

#code for fit with subclonal vaf in exsubcl:
    # excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
    # invf <- 1/exsubcl - 1/higheraf # subtract to use fit without intercept
    # maxi <- length(invf)
    # labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/1.5), maxi)
    # model <- lm(excum~invf+0)
    # sfit <- summary(model)
    # print(sfit)
    # coeffs <- coefficients(model)
    # beta <- unname(coeffs[2])
    # int <- unname(coeffs[1])
    # dr2 <- data.frame(r=sfit$r.squared, intercept=int, slope=beta, subcl=length(exsubcl), all=length(exsubcl_nohigh))
    # write.table(dr2, file=r2, sep="\t", quote=F)
    # pdf(fit)
    # plot(invf, excum, cex=1.5, xaxt="n", xlab='1/f', ylab="Cumulative n. of muts M(f)", main=paste0("R2=", round(sfit$r.squared, digits=3)))
    # oi <- invf[order(invf)]
    # oex <- exsubcl[order(-exsubcl)]
    # axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
    # abline(model, col="red")

excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]>=exsubcl[1:length(exsubcl)]))
excum_p <- excum
exsubcl_p <- exsubcl


load(lm_fit)

excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]>=exsubcl[1:length(exsubcl)]))
data <- data.frame(n=c(excum_p, excum), VAF=c(exsubcl_p, exsubcl), type=c(rep('CRC1599PR', length(excum_p)), rep('CRC1599LM', length(excum))))
x_breaks<-guess_ticks(data$f, fixed_min=0.12, fixed_max=0.24)
y_breaks<-guess_ticks(data$n, fixed_max=68)

data$type <- factor(data$type, levels= c('CRC1599PR', 'CRC1599LM'))

p <- ggplot(data=data, aes(x=VAF, y=n, color=type))+geom_point(size=1)+
  xlab('VAF')+ ylab("Cumulative n. of muts < VAF")+scale_color_manual(values=c('#ffcc33', '#ff9900'))+
  unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand = c(0, 0))+
  theme(legend.position="none",  legend.spacing.y = unit(0.15, "mm")) 
print(p)
graphics.off()
ggsave(outplot, plot=p, width=60, height=60, units="mm")


save.image(paste0(outplot, '.Rdata'))