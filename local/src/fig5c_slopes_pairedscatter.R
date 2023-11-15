data_f  <- snakemake@input[['data']]
log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))
library(ggplot2)
library(ggpubr)
load(theme)

data <- read.table(data_f, sep="\t", header=TRUE)

data$lmodel <- substr(rownames(data), 0, 10)
data$smodel <- substr(rownames(data), 0, 7)
data$mp <- substr(rownames(data), 8, 10)
data$mp <- factor(data$mp, levels=c('PRX', 'LMX'))

fit_r2 <- data
fit_keep <- fit_r2[fit_r2$r > 0.90 & fit_r2$subcl > 10,]
n <- nrow(fit_keep)

pairs <- as.data.frame(table(fit_keep$smodel))
with_pair <- fit_keep[fit_keep$smodel %in% pairs[pairs$Freq == 2,'Var1'],]
fit_r2 <- with_pair
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  #error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  error <- vec_sd
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}


LEVEL <- 0.99
ic_clones <- sapply(unique(fit_r2$mp), function(x) { confidence_interval(fit_r2[fit_r2$mp==x,'intercept'], LEVEL) })
colnames(ic_clones) <- unique(fit_r2$mp)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)
y_breaks<-guess_ticks(fit_r2$intercept,fixed_max=64)
pdata$model <- factor(pdata$model, levels=c('PRX', 'LMX'))
fit_r2$model <- factor(fit_r2$mp, levels=c('PRX', 'LMX'))
pd <- position_dodge(width=0.2)
y_breaks<-guess_ticks(fit_r2$intercept,fixed_max=68)
print(pdata)
p <- ggplot(data=fit_r2, aes(x=model, y=intercept, color=mp, group=smodel)) +
  geom_jitter(data=fit_r2, aes(x=mp, y=intercept, color=mp), size=2, shape=18, position=pd)+
  geom_line(data=fit_r2, aes(group=smodel), position=pd, color="lightgrey", linetype = "dashed",)+
  geom_point(data=pdata, aes(x=model, y=mean, group=NULL, color=NULL), stat="identity", shape=1, size=2) +
  #geom_segment(data=pdata, aes(y=lower, yend=upper, x=model, xend=model, group=NULL, color=NULL), size=0.6) +
  geom_errorbar(data=pdata, aes(x=model, y=lower, ymin=lower, ymax=upper,group=NULL, color=NULL), width=.1, size=.2)+
  scale_color_manual(values=c('#adacac', '#595959'))+
  #scale_y_continuous(breaks = seq(0, 70, by = 10))+
  ylab('Z')+xlab('PDXs')+unmute_theme+theme(legend.position="none")+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))
  #scale_y_continuous(breaks=y_breaks)

  #geom_signif(data=fit_r2, mapping=aes(x=mp, y=intercept), 
  #            comparisons = list(c("LMX", "PRX")), test="t.test", test.args=list(alternative = "greater", paired=TRUE))
pdf('fig_5c_slopes_pairedscatter.pdf')
print(p)
graphics.off()

ggsave(outplot, plot=p, width=60, height=60, units="mm")


pp <- p + theme(legend.position= "none")
ggsave(paste0('nolegend_', outplot), plot=pp, width=60, height=60, units="mm")


met <- fit_r2[fit_r2$mp == "LMX",]
pri <- fit_r2[fit_r2$mp == "PRX",]
if (!all(met$smodel==pri$smodel)) {
  stop('llama! Qualquadra non cosa in pri-met pairs')
}
wt <- wilcox.test(met$intercept, pri$intercept, alternative="greater", paired=TRUE)

sink(log_f)
print('num samples > thr R2 subcl')
print(n)
print('num pri')
length(pri$intercept)
median(pri$intercept)
print('num met')
length(met$intercept)
median(met$intercept)
print(wt)
table(met$intercept > pri$intercept)
sink()

save.image(paste0(outplot, '.Rdata'))

