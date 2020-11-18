#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
inputf <- args[1]
outputf1 <- args[2]
outputf2 <- args[3]


library(ggplot2)
library(dplyr)
#setwd('/mnt/trcanmed/snaketree/prj/AF_spectra/dataset/CRC1307_platypus_nobin')

data <- read.table(inputf, header=TRUE, sep="\t")
vdata <- data[data$class==1,]
vdata$clone <- as.factor(vdata$clone)

gd <- vdata %>%
group_by(clone) %>%
summarise(
MR_EDU = mean(MR_EDU),
MR_conte = mean(MR_conte)
)

ggplot(vdata, aes(x = clone, y = MR_EDU, color=clone, fill=clone)) + geom_point(alpha = .7, size=1, shape=24) + geom_point(data = gd, size = 4)+theme_bw()
ggsave(outputf1)


# TODO FIXME
#vdata <- data[grepl('M', data$clone, fixed=TRUE),]

gd$condition <- 'invivo'
vdata$condition <- ifelse(vdata$class==1,'invitro','invivo')
vdata$condition <- as.factor(vdata$condition)

gd <- vdata %>%
  group_by(condition) %>%
  summarise(
    MR_EDU = mean(MR_EDU),
    MR_conte = mean(MR_conte)
  )

ggplot(vdata, aes(x = condition, y = MR_EDU, color=condition, fill=condition)) + geom_point(alpha = .7, size=1, shape=24) + geom_point(data = gd, size = 4)+theme_bw()

ggsave(outputf2)

