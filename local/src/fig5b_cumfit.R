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

invf_p <- invf
excum_p <- excum
model_p <- model
coeffs_p <- coefficients(model_p)
beta_p <- unname(coeffs_p[1]) # it's 1 because it's a fit without intercept

load(lm_fit)
coeffs <- coefficients(model)
beta <- unname(coeffs[1])

data <- data.frame(n=c(excum_p, excum), invf=c(invf_p, invf), type=c(rep('PRX', length(excum_p)), rep('LMX', length(excum))))
#x_breaks<-guess_ticks(data$invf)
y_breaks<-guess_ticks(data$n, fixed_max=68)

data$type <- factor(data$type, levels= c('PRX', 'LMX'))

# labels on x will probably be manually added at:
# 1/0.25 1/0.2 1/0.18 1/0.15 1/0.12
# labels will need to be not the numbers but the text of the ratios.
#save.image('cum.Rdata')
x_breaks <- c(1/0.25, 1/0.2, 1/0.15, 1/0.12) - 1/higheraf
x_labels <- c('1/0.25','1/0.2','1/0.15','1/0.12')
p <- ggplot(data=data, aes(x=invf, y=n, color=type))+geom_point(size=1)+#stat_smooth(method = "lm", formula="excum~invf+0", se=FALSE, show.legend = FALSE, size=0.5)+ 
  geom_abline(slope=beta, intercept=0, color='#ff9900', size=0.5)+
  geom_abline(slope=beta_p, intercept=0, color='#ffcc33', size=0.5)+
  xlab('1/f')+ ylab("Cumulative n. of muts M(f)")+scale_color_manual(values=c('#ffcc33', '#ff9900'))+
  unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, labels=x_labels, limits=c(min(x_breaks),max(x_breaks)), expand = c(0, 0))+
  theme(legend.position="none",  legend.spacing.y = unit(0.15, "mm")) 
pdf('fig_5b_cumfit.pdf')
print(p)
graphics.off()
ggsave(outplot, plot=p, width=60, height=60, units="mm")


save.image(paste0(outplot, '.Rdata'))