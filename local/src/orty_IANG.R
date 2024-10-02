
library(ggplot2)
load('/scratch/trcanmed/AF_spectra/dataset_fIANG/fig_1b_MR.svg.Rdata')

our <- our[grepl('PR', our$sample),]
pdata <- pdata[grepl('PR', pdata$model),]
our$model <- droplevels(our$model)
pdata$xmodel <- c(4,2,3,1)
y_breaks <- guess_ticks(our$MR)
pal <- pal[grepl('PR', names(pal))]
p <- ggplot() + 
  geom_point(data=our, aes(x=model, y=MR, color=model_clone), stat="identity", size=2, shape=18, position=position_dodge(0.7))+
  geom_segment(data=pdata, aes(x=xmodel-0.2, yend=mean,y=mean,  xend=xmodel+0.2),size=.3) +
  geom_errorbar(data=pdata, aes(ymin=lower, ymax=upper, x=model), size=0.3, width=0.3)+ylab('MR [SNV/(Gbp*division)]')+xlab('PDTs')+
  scale_color_manual(values=pal)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, 1.5),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   
#print(p)
#graphics.off()
#print(our)
ggsave('focus_PRI.svg', plot=p, width=89, height=89, units="mm")