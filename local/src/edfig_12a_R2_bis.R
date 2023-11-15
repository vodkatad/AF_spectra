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

data <- data[data$laf != 0.025,]
data <- data[data$haf != 0.050,]

y_breaks<-c(0,0.25,0.5,0.75,1)
#x_breaks<-as.factor(unique(data$laf))

y_labels <- y_breaks
#y_labels[1] <- '-18.0'
r2_plot <- ggplot(data=data, aes(x=as.factor(laf), y=r)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=0.1, alpha=0.7, height=0) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "red")+
  facet_grid(~haf) +
  xlab("Lower AF")+ylab("R")+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0), labels=y_labels)+
  theme_bw()+
  theme(
	text = element_text(size = textSize, family='sans'),
	axis.title = element_text(size = largerSize),
	axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
	axis.text.y = element_text(size = textSize, color="black"),
	plot.title = element_text(size = largerSize, hjust = 0.5),
	legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
	axis.line = element_line(colour = "black"),
	axis.ticks = element_line(color = "black")
)
  
get_delta_all <- function(laf, haf, data) {
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
  #slopeP <- mean(pri$intercept)
  #slopeM <- mean(met$intercept)
  #meanD <- mean(met$intercept) - mean(pri$intercept)
  if (nrow(met) > 0) {
    return(met$intercept - pri$intercept)
  } else {
    return(NA)
  }
}

dd <- as.data.frame(table(data$laf, data$haf))
dd <- dd[dd$Freq != 0,]

deltas <- mapply(get_delta_all, as.numeric(as.character(dd$Var1)), as.numeric(as.character(dd$Var2)), MoreArgs=list(data))

delta_df <- data.frame(delta=c(), haf=c(), laf=c())
laf <- as.numeric(as.character(dd$Var1))
haf <- as.numeric(as.character(dd$Var2))

for (i in seq(1,length(deltas))) {
  myd <- deltas[[i]]
  add <- data.frame(delta=myd, haf=rep(haf[i], length(myd)), laf=rep(laf[i], length(myd)))
  delta_df <- rbind(delta_df, add)
}

## asterisks
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

ma <- mapply(get_info_aft, as.numeric(as.character(dd$Var1)), as.numeric(as.character(dd$Var2)), MoreArgs=list(data))
topl <- as.data.frame(t(ma))
colnames(topl) <- c('n_models', 'meanR2', 'slopeP', 'slopeM', 'deltasl', 'tpval')
topl$name <- paste(as.numeric(as.character(dd$Var1)), as.numeric(as.character(dd$Var2)), sep='-')
topl$laf <- as.numeric(as.character(dd$Var1))
topl$haf <- as.numeric(as.character(dd$Var2))

topl$sign <- ifelse(topl$tpval < 0.05,'Yes', 'No')
topl$asterisk <- ifelse(topl$tpval <= 0.0001, '****', ifelse(topl$tpval <= 0.001, '***', 
                  ifelse(topl$tpval <= 0.01, '**', ifelse(topl$tpval <= 0.05, '*', ' ')) ))       



topl$y <- rep(44, nrow(topl))

y_breaks <- guess_ticks(delta_df$delta, fixed_min=-15, fixed_max=45)
#y_breaks <-  guess_ticks(delta_df$delta, fixed_min=-20, fixed_max=60)
y_labels <- y_breaks
y_labels[1] <- 0.25 # trick to have y axis aligned
deltas_plot <- ggplot(data=delta_df, aes(x=as.factor(laf), y=delta)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=0.1, alpha=0.7, height=0) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "red")+
  facet_grid(~haf) +
  xlab("Lower AF")+ylab("Zm-Zp")+
  scale_y_continuous(breaks=y_breaks,limits=c(min(y_breaks), max(y_breaks)),expand = c(0, 0), labels=y_labels)+
  geom_text(data=topl, aes(y=y, x=as.factor(laf), label=asterisk))+
  theme_bw()+
  theme(
    text = element_text(size = textSize, family='sans'),
    axis.title = element_text(size = largerSize),
    axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
    axis.text.y = element_text(size = textSize, color="black"),
    plot.title = element_text(size = largerSize, hjust = 0.5),
    legend.title = element_text(size=largerSize, hjust = 0.5),
    legend.text = element_text(size=textSize),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(color = "black")
  )

p <- ggarrange(deltas_plot, r2_plot, 
          labels = c("", ""),
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