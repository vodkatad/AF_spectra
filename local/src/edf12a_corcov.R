MR_f  <- snakemake@input[['MR']]
MRcov_f  <- snakemake@input[['MRcov']]
colors <- snakemake@input[['palette']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

name1 <- "MR1x"
name2 <- "MR2x"

data1 <- read.table(MR_f, sep="\t", stringsAsFactors = FALSE)
data2 <- read.table(MRcov_f, sep="\t", stringsAsFactors = FALSE)
colnames(data1) <- c('sample', name1)
colnames(data2) <- c('sample', name2)

m <- merge(data1, data2, by="sample")

m[, name1] <- m[, name1] / 0.000000001
m[, name2] <- m[, name2] / 0.000000001
m$model <- sapply(m$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
m$clone <- sapply(m$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
m$clone2 <- sapply(m$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
m$model_clone <- paste0(m$model, "_", m$clone)
m$time <- sapply(m$sample, function(x) {y<-strsplit(x, '-')[[1]][3]; return(y[1])})

models <- c('CRC1307',  'CRC1502')
sink(log_f)
print('Spearman on single clones estimates')
print(models)
sink()
for (mi in models) {
  sink(log_f, append=TRUE)
  onem <- m[m$model == mi,]
  print(mi)
  print(cor.test(onem[, name1], onem[, name2], method="spearman"))
  sink()
}

#ggplot(data=m, aes(x=platy, y=noplaty, color=model_clone))+geom_point()+
#  unmute_theme+scale_color_manual(values=pal)

death_conversion_dpi96 <- 96/72
print(m)
y_breaks <- guess_ticks(m$MR2x)
x_breaks <- guess_ticks(m$MR1x)
p <- ggplot(data=m, aes_string(x=name1, y=name2, color='model_clone'))+geom_point(size=1, shape=16, alpha=0.8)+
  geom_abline(slope=1, intercept=0)+
  unmute_theme+scale_color_manual(values=pal)+xlab('MR on >1x regions')+ylab('MR on >20x regions')+
  scale_y_continuous(breaks=y_breaks, limits=c(0,max(y_breaks)), expand = c(0, 0))+
  scale_x_continuous(breaks=x_breaks, limits=c(0,max(x_breaks)), expand = c(0, 0))+
  theme(legend.position="none")
ggsave(outplot, p, height=60*death_conversion_dpi96, width=60*death_conversion_dpi96, units="mm")

prepare_plot <- function(our) {
  confidence_interval <- function(vector, interval) {
    # Standard deviation of sample
    vec_sd <- sd(vector)
    # Sample size
    n <- length(vector)
    # Mean of sample
    vec_mean <- mean(vector)
    # Error according to t distribution
    error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
    # Confidence interval as a vector
    result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
    return(result)
  }
  
  LEVEL <- 0.99
  ic_clones <- sapply(unique(our$model), function(x) { confidence_interval(our[our$model==x,'MR'], LEVEL) })
  colnames(ic_clones) <- unique(our$model)
  pdata <- as.data.frame(t(ic_clones))
  pdata$model <- rownames(pdata)
  return(pdata)
}

our1 <- m
our1[,name2] <- NULL
colnames(our1)[colnames(our1) == name1] <- 'MR'

our2 <- m
our2[,name1] <- NULL
colnames(our2)[colnames(our2) == name2] <- 'MR'

pdata1 <- prepare_plot(our1)
pdata2 <- prepare_plot(our2)
pdata1$class <- name1
pdata2$class <- name2

our1$class <- name1
our2$class <- name2

our <- rbind(our1, our2)
pdata <- rbind(pdata1, pdata2)

sink(log_f, append=TRUE)
print(pdata1$mean / pdata2$mean)
print('total clones')
print(nrow(m))
print(table(m$model))
sink()

save.image(paste0(outplot, '.Rdata'))