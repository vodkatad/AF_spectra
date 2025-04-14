load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/fig_3a_MR.svg.Rdata')

wilcox.test(formula=as.formula('MR_edu~time'), data=our[our$model == "CRC0282PR", ])

wilcox.test(formula=as.formula('MR_edu~time'), data=our[our$model == "CRC1502LM", ])

wilc <- function(model, data) {
  w <- wilcox.test(formula=as.formula('MR_edu~time'), data=data[data$model == model,])
  return(w$p.value)
}


ww <- as.data.frame(sapply(unique(our$model), wilc, our))
colnames(ww) <- 'p'
ww$padj <- p.adjust(ww$p)


wilc_mc <- function(model, data) {
  smodel <- substr(model, 0, 7)
  d <- data[grepl(smodel, data$model),]
  d <- d[grepl(model, d$model_clone),]
  if (length(unique(d$time)) == 2) {
    w <- wilcox.test(formula=as.formula('MR_edu~time'), data=d)
    return(w$p.value)
  } else {
    return(NA)
  }
}

pval <- as.data.frame(sapply(unique(our$model_clone), wilc_mc, our))
colnames(pval) <- 'p'
pval <- pval[!is.na(pval$p),, drop=F]

pval$padj <- p.adjust(pval$p)

#Possiamo provare un approccio Bayesiano in cui, per ogni misura a T2, stimiamo se è più probabile che appartenga alla distribuzione di misure a T1 oppure a una distribuzione indipendente?


shap <- function(model, data) {
  sh <- shapiro.test(data[data$model == model, 'MR_edu'])
  return(sh$p.value)
}

sha <- as.data.frame(sapply(unique(our$model), shap, our))
colnames(sha) <- 'p'
sha$padj <- p.adjust(sha$p)

##

shap_time <- function(model, data, time) {
  sh <- shapiro.test(data[data$model == model & data$time == time, 'MR_edu'])
  return(sh$p.value)
}

sha_1 <- as.data.frame(sapply(unique(our$model), shap_time, our, 1))
colnames(sha_1) <- 'p'
sha_1$padj <- p.adjust(sha_1$p)

sha_2 <- as.data.frame(sapply(unique(our$model), shap_time, our, 2))
colnames(sha_2) <- 'p'
sha_2$padj <- p.adjust(sha_2$p)


# anova
aa <- aov(MR_edu ~ time * model,
            data = our,
            var.equal = TRUE # assuming equal variances
)
summary(aa)


library(multcomp)

# Tukey HSD test:

post_test <- glht(aa)

#                  , linfct = mcp(model = "Tukey"))

summary(post_test)
#  shapiro.test(aa$residuals)  mh

# library(rstatix)
# res.aov <- anova_test(data = as_tibble(our), dv = MR_edu, wid = model, within = time)
# get_anova_table(res.aov)

# posttest

# kw


### t test single clones
wilc_mc <- function(model, data) {
  smodel <- substr(model, 0, 7)
  d <- data[grepl(smodel, data$model),]
  d <- d[grepl(model, d$model_clone),]
  if (nrow(d)>4) {
    w <- t.test(formula=as.formula('MR_edu~time'), data=d)
    return(w$p.value)
  } else {
    return(NA)
  }
}

pval <- as.data.frame(sapply(unique(our$model_clone), wilc_mc, our))
colnames(pval) <- 'p'
pval <- pval[!is.na(pval$p),, drop=F]

pval$padj <- p.adjust(pval$p)

#
library(BayesFactor)
bay <- function(model, data) {
  w <- ttestBF(formula=as.formula('MR_edu~time'), data=data[data$model == model,])
  print(w)
}


bb <- as.data.frame(sapply(unique(our$model), bay, our))
