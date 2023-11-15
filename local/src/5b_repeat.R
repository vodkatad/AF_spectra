load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_5b_cumfit.svg.Rdata')

# Cumulative number <= VAF, vs VAF
load(pr_fit)

excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]>=exsubcl[1:length(exsubcl)]))

invf_p <- invf
exsubcl_p <- exsubcl
excum_p <- excum
model_p <- model
coeffs_p <- coefficients(model_p)
beta_p <- unname(coeffs_p[1]) # it's 1 because it's a fit without intercept

load(lm_fit)
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]>=exsubcl[1:length(exsubcl)]))

coeffs <- coefficients(model)
beta <- unname(coeffs[1])

data <- data.frame(n=c(excum_p, excum), f=c(exsubcl_p, exsubcl), type=c(rep('CRC1599PR', length(excum_p)), rep('CRC1599LM', length(excum))))
x_breaks<-guess_ticks(data$f, fixed_min=0.12, fixed_max=0.24)
y_breaks<-guess_ticks(data$n, fixed_max=68)

data$type <- factor(data$type, levels= c('CRC1599PR', 'CRC1599LM'))

# labels on x will probably be manually added at:
# 1/0.25 1/0.2 1/0.18 1/0.15 1/0.12
# labels will need to be not the numbers but the text of the ratios.
#save.image('cum.Rdata')
#x_breaks <- c(1/0.25, 1/0.2, 1/0.15, 1/0.12) - 1/higheraf
#x_labels <- c('1/0.25','1/0.2','1/0.15','1/0.12')
p <- ggplot(data=data, aes(x=f, y=n, color=type))+geom_point(size=1)+
 xlab('VAF')+ ylab("Cumulative n. of muts > VAF")+scale_color_manual(values=c('#ffcc33', '#ff9900'))+
  #unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand = c(0, 0))+
  theme_bw()
