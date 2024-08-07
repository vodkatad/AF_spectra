R2_slopes_f  <- snakemake@input[['data']]
log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
load(theme)

library(stringr)
library(reshape)
library(ggplot2)
library(ggpubr)
save.image(paste0(outplot, '.Rdata'))

data <- read.table(R2_slopes_f, sep="\t", header=TRUE)
data$lmodel <- substr(rownames(data), 0, 10)
data$smodel <- substr(rownames(data), 0, 7)
data$mp <- substr(rownames(data), 8, 10)
data$mp <- factor(data$mp, levels=c('PRX', 'LMX'))

get_laf <- function(s) {
  as.numeric(str_match(s, "fit.(\\d+.\\d+)")[1,2])
}

get_haf <- function(s) {
  as.numeric(str_match(s, "(\\d+.\\d+).r2")[1,2])
}

data$laf <- sapply(rownames(data), get_laf)
data$haf <- sapply(rownames(data), get_haf)
#data <- data[data$laf != 0.025,]


r2_plot <- ggplot(data=data, aes(x=as.factor(laf), y=r)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=0.1, alpha=0.7) +
  facet_grid(~haf) +
  xlab("Lower AF")+ylab("R²")+
  unmute_theme


get_info_aft <- function(laf, haf, data) {
  sub <- data[data$laf == laf & data$haf == haf,]
  fit_keep <- sub[sub$r > 0.90 & sub$subcl > 10,]
  pairs <- as.data.frame(table(fit_keep$smodel))
  with_pair <- fit_keep[fit_keep$smodel %in% pairs[pairs$Freq == 2,'Var1'],]
  met <- with_pair[with_pair$mp == "LMX",]
  pri <- with_pair[with_pair$mp == "PRX",]
  if (!all(met$smodel==pri$smodel)) {
    stop('llama! Qualquadra non cosa in pri-met pairs')
  }
  tp <- NA
  if (nrow(met) > 1 & nrow(pri) > 1) {
    ti <- wilcox.test(met$intercept, pri$intercept, alternative="greater", paired=TRUE)
    tp <- ti$p.value
  }
  slopeP <- mean(pri$intercept)
  slopeM <- mean(met$intercept)
  meanD <- mean(met$intercept) - mean(pri$intercept)
  return(c(length(unique(with_pair$smodel)), mean(with_pair$r), slopeP, slopeM, meanD, tp))
}

dd <- as.data.frame(table(data$laf, data$haf))
dd <- dd[dd$Freq != 0,]

ma <- mapply(get_info_aft, as.numeric(as.character(dd$Var1)), as.numeric(as.character(dd$Var2)), MoreArgs=list(data))

topl <- as.data.frame(t(ma))
colnames(topl) <- c('n_models', 'meanR2', 'slopeP', 'slopeM', 'deltasl', 'tpval')
topl$name <- paste(as.numeric(as.character(dd$Var1)), as.numeric(as.character(dd$Var2)), sep='-')
topl$laf <- as.numeric(as.character(dd$Var1))
topl$haf <- as.numeric(as.character(dd$Var2))

topl$sign <- ifelse(topl$tpval < 0.05,'Yes', 'No')
delta_plot <- ggplot(data=topl, aes(x=as.factor(laf), y=deltasl, fill=sign)) + 
  geom_col() +
  facet_grid(~haf) +
  scale_fill_manual(values=c('darkgreen', 'darkgoldenrod'))+
  ylab('Slope LMX - Slope PRX') +
  xlab('Lower VAF') +
  geom_text(aes(label=n_models), position=position_dodge(width=0.9), vjust=-0.15, size=2) +
  theme(legend.position= "none")+unmute_theme
  

p <- ggarrange(r2_plot, delta_plot, 
          labels = c("a", ""),
          ncol = 1, nrow = 2)

ggsave(outplot, plot=p, width=183, height=90, units="mm")

#ggsave('Redo_mixedPriMets.pdf', width=11.68, height=8.26, units="in")

#topl$sign <- topl$tpval < 0.05
#ggplot(data=topl, aes(x=as.factor(laf), y=deltasl, fill=sign)) + 
#  geom_col() +
#  theme_bw() +
#  facet_grid(~haf) +
#  theme(text=element_text(size=18))


# ggplot(data=topl, aes(x=as.factor(laf), y=n_models, fill=sign)) + 
#   geom_col() +
#   theme_bw() +
#   facet_grid(~haf) +
#   theme(text=element_text(size=18))

# wideR2 <- cast(data=topl, formula="laf ~ haf", value="meanR2") 
# widesl <- cast(data=topl, formula="laf ~ haf", value="deltasl")
# rownames(wideR2) <- wideR2$laf
# wideR2$laf <- NULL
# rownames(widesl) <- widesl$laf
# widesl$laf <- NULL

# topl$d <- round(topl$deltasl, digits=3)
# plot(data=topl, aes(x=as.factor(laf), y=as.factor(haf), fill=meanR2))+geom_tile()+theme_bw()+
# geom_text(data=topl,aes(x=as.factor(laf), y=as.factor(haf), label=d))+
# scale_fill_continuous(type = "viridis")+theme(text=element_text(size=18))

# ggsave(outplot, plot=p, width=60, height=60, units="mm")

save.image(paste0(outplot, '.Rdata'))