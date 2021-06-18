library(ggplot2)
th <- function() {
  textSize <- 1.5
  current_theme <-
    theme_bw() +
    theme(
      strip.text = element_text(size = rel(textSize)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title = element_text(size = rel(1.8)),
      axis.text.x = element_text(size=rel(1.7)),
      axis.text.y = element_text(angle = 0,
                                 size = rel(1.7)),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.length.y.left = unit(3,'mm'),
      
      legend.position = "top",
      legend.justification = "right",
      #legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(size = rel(textSize), face = "bold"),
      legend.text = element_text(size = rel(1.2)),
      legend.background = element_rect(size=0.5, linetype="solid", color="black"),
      plot.title = element_text(
        face = "bold",
        size = rel(2),
        hjust = 0.5
      ),
      panel.border = element_blank(),
      plot.caption = element_text(size=rel(1))
    )
  current_theme
}
current_theme <- th()

setwd('/scratch/trcanmed/AF_spectra/dataset/CRC1307_clones_all/tree')
load('tree_bulk_vitro.Rdata')
library(ape)
di <- cophenetic(NexusTree)

names <-rowNames(di)

# TODO automagicate (from conf?)
names_tree <- list(
  'CRC1307-02-1-A'= 'CRC1307-02-0',
  'CRC1307-02-1-B'= 'CRC1307-02-0',
  'CRC1307-02-1-E'= 'CRC1307-02-0',
  'CRC1307-08-1-B'= 'CRC1307-08-0',
  'CRC1307-08-1-D'= 'CRC1307-08-0',
  'CRC1307-08-1-E'= 'CRC1307-08-0',
  'CRC1307-09-1-B'= 'CRC1307-09-0',
  'CRC1307-09-1-C'= 'CRC1307-09-0',
  'CRC1307-09-1-E'= 'CRC1307-09-0',
  'CRC1307-09E-2-3'= 'CRC1307-09-1-E',
  'CRC1307-09E-2-4'= 'CRC1307-09-1-E',
  'CRC1307-09E-2-5'= 'CRC1307-09-1-E'
) 

distances <- sapply(names(names_tree), function(x) { di[rownames(di)==x, colnames(di)==names_tree[[x]]]} )

tree_di <- as.data.frame(distances)

mr <- read.table('/scratch/trcanmed/AF_spectra/dataset/MR_edu_SNV', sep="\t", header=FALSE, row.names = 1)
colnames(mr) <- c('MR')

m_mr <- merge(tree_di, mr, by="row.names")

p1 <- ggplot(data=m_mr, aes(x=distances, y=MR))+geom_point()+current_theme+stat_smooth(method="lm")

g <- read.table('/scratch/trcanmed/AF_spectra/dataset/vitro_gained_SNV', sep="\t", header=FALSE, row.names = 1)
colnames(g) <- c('gained')
m_g <- merge(tree_di, g, by="row.names")
p2 <- ggplot(data=m_g, aes(x=distances, y=gained))+geom_point()+current_theme+stat_smooth(method="lm")




gn <- read.table('/scratch/trcanmed/AF_spectra/dataset/vitro_gained_norm_SNV', sep="\t", header=FALSE, row.names = 1)
colnames(gn) <- c('gainednorm')
m_gn <- merge(tree_di, gn, by="row.names")
p3 <- ggplot(data=m_gn, aes(x=distances, y=gainednorm))+geom_point()+current_theme+stat_smooth(method="lm")

p2

summary(lm(data=m_g, formula="gained~distances"))

cor.test(m_g$gained,m_g$distances)

