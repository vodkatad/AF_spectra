#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = T)
dir <- args[1]

wd <- getwd()
setwd(dir)

objsf <- as.list(c('CRC0282-01-0.AF.rda','CRC0282-05-0.AF.rda','CRC0282-07-0.AF.rda',
                'CRC0327-02-0.AF.rda','CRC0327-04-0.AF.rda','CRC0327-08-0.AF.rda',
                'CRC0441-01-0.AF.rda','CRC0441-03-0.AF.rda','CRC0441-10-0.AF.rda',
                'CRC1078-02-0.AF.rda','CRC1078-07-0.AF.rda','CRC1078-09-0.AF.rda',
                'CRC1307-02-0.AF.rda','CRC1307-08-0.AF.rda','CRC1307-09-0.AF.rda',
                'CRC1502-03-0.AF.rda','CRC1502-08-0.AF.rda','CRC1502-09-0.AF.rda',
                'CRC1502-10-0.AF.rda','CRC1599LM-01-0.AF.rda',
                'CRC1599LM-03-0.AF.rda','CRC1599LM-07-0.AF.rda',
                'CRC1599PR-01-0.AF.rda',
                'CRC1599PR-10-0.AF.rda'))

tot <- length(objsf)

i <- 1
k <- 1
while (i <= tot) {
  is <- seq(i, i+4)
  mylist <- list()
  for (j in seq(1, length(is))) {
    if (is[j] <= tot) {
      mylist <- append(mylist, objsf[[is[j]]])
    }
  }
  ggl <- lapply(mylist, readRDS)
  m <- ggarrange(plotlist=ggl, ncol=1, nrow=5)
  setwd(wd)
  pdf(paste0('temp_S1_soccia_', k, '.pdf'), width=8.27, height=11.69)
  print(m+ theme(plot.margin = margin(t = 58.5, r = 15, b = 58.5, l = 15, unit = "mm")))
  #print(m)
  graphics.off()  
  setwd(dir)
  i <- i + 5
  k <- k + 1
}

#egrassi@qui:/tmp/plots$ pdftk temp_S1_soccia_*.pdf cat output test.pdf


