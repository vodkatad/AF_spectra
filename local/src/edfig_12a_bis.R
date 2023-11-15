load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/edfig_12a_R2.svg.Rdata')

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

deltas <- mapply(get_delta_all, as.numeric(as.character(dd$Var1)), as.numeric(as.character(dd$Var2)), MoreArgs=list(data))

delta_df <- data.frame(delta=c(), haf=c(), laf=c())
laf <- as.numeric(as.character(dd$Var1))
haf <- as.numeric(as.character(dd$Var2))

for (i in seq(1,length(deltas))) {
  myd <- deltas[[i]]
  add <- data.frame(delta=myd, haf=rep(haf[i], length(myd)), laf=rep(laf[i], length(myd)))
  delta_df <- rbind(delta_df, add)
}

y_breaks <- guess_ticks(delta_df$delta, fixed_min=-14, fixed_max=43)
deltas_plot <- ggplot(data=delta_df, aes(x=as.factor(laf), y=delta)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=0.1, alpha=0.7, height=0) +
  facet_grid(~haf) +
  xlab("Lower AF")+ylab("Zm-Zp")+
  scale_y_continuous(breaks=y_breaks,limits=c(min(y_breaks), max(y_breaks)),expand = c(0, 0))+
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