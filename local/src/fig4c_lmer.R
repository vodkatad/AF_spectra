MR_f  <- snakemake@input[['mr']]
CN_f  <- snakemake@input[['cn']]
colors <- snakemake@input[['palette']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
data_f <- snakemake@output[['avgdata']]

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
#library(lme4)
library(nlme)
library(ggeffects)
load(theme)
guess_ticks_underzero <- function(values, nticks=5, fixed_max=NULL) {
  vmax <- max(values)
  if (is.null(fixed_max)) { 
    round_max <- vmax
  } else {
    round_max <- fixed_max
  }
  v_min<-min(values)
  my_breaks <- seq(v_min, round_max, length.out=nticks)
  return(my_breaks)
}
palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model
names(pal) <- paste0(names(pal), ifelse(!grepl('\\d$', names(pal)), '', ifelse(names(pal)=="CRC0282", 'PR', 'LM')))

mr <- read.table(MR_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(mr) <- paste0(colnames(mr), "_mr")
cn <- read.table(CN_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(cn) <- paste0(colnames(cn), "_cn")

m <- merge(mr, cn, by="row.names")
m$model <- m$Row.names
m$mean_mr<-m$mean_mr*1000000000
y_breaks <- guess_ticks(m$mean_cn+6,fixed_max=40)
x_breaks<-guess_ticks(m$mean_mr)

m$patient <- ifelse(m$model %in% c('CRC1307LM', 'CRC1078LM'), 'CRC1307all', substr(m$model, 0, 7))


# checked that it's the same in terms of significance without *10**9

#https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode
#lmer(value~status+(1|experiment)))
#m1 <- lme(value~status,random=~1|experiment,data=mydata)
#mixedlm <- lmer(mean_cn ~ mean_mr + (1|patient), data = m) 

mixedlm <- lme(mean_cn~mean_mr, random=~1|patient, data = m) 
sink(log_f)
summary(mixedlm)
sink()

#https://ourcodingclub.github.io/tutorials/mixed-models/
pred.mm <- ggpredict(mixedlm, terms = c("mean_mr"))  # this gives overall predictions for the model

# Plot the predictions 

p <- ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = m,                      # adding the raw data (scaled values)
               aes(x = mean_mr, y = mean_cn, colour = model)) + 
    unmute_theme+theme(legend.position="none") +
    scale_color_manual(values=pal)+xlab('MR [SNVs/(Gbp*division)]')+ylab('MEDICC2 events')+
    scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
    scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)),expand=c(0,0))


print(p)
ggsave(outplot, plot=p, width=89, height=89, units="mm")

sink(log_f, append=TRUE)
print('n models')
print(nrow(m))
sink()
save.image(paste0(outplot, '.Rdata'))

m <- m[m$model != "CRC0282PR",]
y_breaks <- guess_ticks(m$mean_cn+7)


mixedlm <- lme(mean_cn~mean_mr, random=~1|patient, data = m) 
sink(log_f, append=TRUE)
summary(mixedlm)
sink()

#https://ourcodingclub.github.io/tutorials/mixed-models/
pred.mm <- ggpredict(mixedlm, terms = c("mean_mr"))  # this gives overall predictions for the model

# Plot the predictions 

p <- ggplot(pred.mm) + 
  geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = m,                      # adding the raw data (scaled values)
             aes(x = mean_mr, y = mean_cn, colour = model)) + 
  unmute_theme+theme(legend.position="none") +
  scale_color_manual(values=pal)+xlab('MR [SNVs/(Gbp*division)]')+ylab('MEDICC2 events')+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)),expand=c(0,0))

ggsave(paste0('noMSI_', outplot), plot=p, width=89, height=56, units="mm")
