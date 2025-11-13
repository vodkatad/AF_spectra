library(ggplot2)
library(ggpubr)
setwd('/scratch/trcanmed/AF_spectra/datasetV2/')

objsf <- as.list(c('CRC0282/mutect_nobin/CRC0282-01-0.AF.rda','CRC0282/mutect_nobin/CRC0282-05-0.AF.rda','CRC0282/mutect_nobin/CRC0282-07-0.AF.rda',
                'CRC0327/mutect_nobin/CRC0327-02-0.AF.rda','CRC0327/mutect_nobin/CRC0327-04-0.AF.rda','CRC0327/mutect_nobin/CRC0327-08-0.AF.rda',
                'CRC0441/mutect_nobin/CRC0441-01-0.AF.rda','CRC0441/mutect_nobin/CRC0441-03-0.AF.rda','CRC0441/mutect_nobin/CRC0441-10-0.AF.rda',
                'CRC1078/mutect_nobin/CRC1078-02-0.AF.rda','CRC1078/mutect_nobin/CRC1078-07-0.AF.rda','CRC1078/mutect_nobin/CRC1078-09-0.AF.rda',
                'CRC1307/mutect_nobin/CRC1307-02-0.AF.rda','CRC1307/mutect_nobin/CRC1307-08-0.AF.rda','CRC1307/mutect_nobin/CRC1307-09-0.AF.rda',
                'CRC1502/mutect_nobin/CRC1502-03-0.AF.rda','CRC1502/mutect_nobin/CRC1502-08-0.AF.rda','CRC1502/mutect_nobin/CRC1502-09-0.AF.rda',
                'CRC1502/mutect_nobin/CRC1502-10-0.AF.rda','CRC1599LM/mutect_nobin/CRC1599LM-01-0.AF.rda',
                'CRC1599LM/mutect_nobin/CRC1599LM-03-0.AF.rda','CRC1599LM/mutect_nobin/CRC1599LM-07-0.AF.rda',
                'CRC1599PR/mutect_nobin/CRC1599PR-01-0.AF.rda',
                'CRC1599PR/mutect_nobin/CRC1599PR-10-0.AF.rda'))

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
  #pdf(paste0('soccia_', k, '.pdf'), width=8.27, height=11.69)
  #print(m+ theme(plot.margin = margin(t = 58.5, r = 15, b = 58.5, l = 15, unit = "mm")))
  #print(m)
  #graphics.off()  
  png(paste0('SupFig.1_page_', k, '.jpg'), width=8.27, height=11.69, res=300, units='in')
  print(m+ theme(plot.margin = margin(t = 58.5, r = 15, b = 58.5, l = 15, unit = "mm")))
  #print(m)
  graphics.off()  
  i <- i + 5
  k <- k + 1
}

#egrassi@qui:/tmp/plots$ pdftk soccia_*.pdf cat output test.pdf


