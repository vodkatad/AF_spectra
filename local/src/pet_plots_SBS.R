library(ggplot2)

textSize <- 15
largerSize <- textSize + 2

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
unmute_theme <- theme(
  text = element_text(size = textSize, family='sans'),
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(color = "black"),
  panel.background = element_blank()
)


slope_1 <- 0.8
slope_8 <- 1.2
int_1 <- 0
int_8 <- -20

x <- seq(0, 100, by=0.5)
y_sbs1 <- x*slope_1+int_1
y_sbs8 <- x*slope_8+int_8

y_sbs8 <- y_sbs8[x >= 50]

pd <- data.frame(x=c(x, x[x>=50]), y=c(y_sbs1, y_sbs8), sbs=c(rep('SBS1', length(y_sbs1)), rep('SBS8', length(y_sbs8))))
ggplot(pd, aes(x=x, y=y, color=sbs))+geom_line(aes(group=sbs), size=1)+
  scale_y_continuous(expand = c(0, 0))+scale_color_manual(values=c('darkgrey', 'coral3'))+
  scale_x_continuous(limits=c(0, 100),expand = c(0, 0),breaks=c(50), labels=c('Cancer initiation'))+unmute_theme+
  theme(axis.ticks=element_blank(), axis.text.y = element_blank(), 
        legend.title = element_blank(), axis.title.x = element_text(hjust = 1, vjust = 0),
        axis.title.y = element_text(hjust = 1, vjust = 0))+
  xlab('time')+ylab('SBS')+
  geom_vline(xintercept=50, linetype = "dashed", linewidth=0.5)

x1 <- 60
x2 <- 60.5

pd2 <- pd[pd$x >= x1 & pd$x <= x2,]

pd2 <- pd[pd$x >= 45 & pd$x <= 65,]

ggplot(pd2, aes(x=x, y=y, color=sbs))+geom_line(aes(group=sbs), size=1)+
  scale_y_continuous(limits=c(40, 60), expand = c(0, 0))+scale_color_manual(values=c('darkgrey', 'coral3'))+
  scale_x_continuous(limits=c(55, 65),expand = c(0, 0), breaks=c(60.0, 60.5), labels=c('T0', 'T1'))+
  unmute_theme+
  theme(axis.ticks=element_blank(), axis.text.y = element_blank(), 
        legend.title = element_blank(), axis.title.x = element_text(hjust = 1, vjust = 0),
        axis.title.y = element_text(hjust = 1, vjust = 0))+  xlab('time')+ylab('SBS')+
  geom_vline(xintercept=60, linetype = "dotted", linewidth=1)+
  geom_vline(xintercept=60.5, linetype = "dotted", linewidth=1)
  
pp <- pd2[pd2$sbs=='SBS8', 'y'] /pd2[pd2$sbs=='SBS1', 'y']

pd3 <- data.frame(t=c('T0', 'T1'), y=pp)

ggplot(pd3, aes(x=t, y=y))+geom_line(size=1, group=0)+geom_point(size=3)+
  scale_y_continuous(limits=c(1.08, 1.09), expand = c(0, 0))+unmute_theme+
  unmute_theme+
  theme(axis.ticks=element_blank(), axis.text.y = element_blank(), 
        legend.title = element_blank(), axis.title.x = element_text(hjust = 1, vjust = 0),
        axis.title.y = element_text(hjust = 1, vjust = 0))+  xlab('time')+ylab('SBS8/SBS1')+
  geom_vline(xintercept=1, linetype = "dotted", linewidth=1)+
  geom_vline(xintercept=2, linetype = "dotted", linewidth=1)


