rdata_f  <- snakemake@input[['rdata']]
my_outplot <- snakemake@output[['plot']]
my_log_f <- snakemake@log[['log']]


load(rdata_f)
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)

#p <- ggplot() + 
#  geom_point(data=our, aes(x=model, y=MR, color=model_clone, group=sample), position=position_dodge(0.8), size=2, shape=18)+
#  geom_segment(data=pdata, aes(x=xmodel-0.2, yend=mean,y=mean,  xend=xmodel+0.2),size=.3) +
#  geom_errorbar(data=pdata, aes(ymin=lower, ymax=upper, x=model), size=0.3, width=0.3)+ylab('MR [SNV/(Gbp*division)]')+xlab('PDTs')+
#  scale_color_manual(values=pal)+
#  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
#  unmute_theme+theme(legend.position="none", axis.text.x = element_blank(), 
#                     axis.ticks.x = element_blank(),
#                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   

#ggsave(outplot, plot=p, width=89, height=89, units="mm")

our$type <- ifelse(grepl('PR', our$sample), 'PRs', 'LMs')
our$type <- factor(our$type, levels=c('PRs', 'LMs'))


pp <- ggplot(data=our, aes(x=type, y=MR))+geom_boxplot(outlier.shape=NA)+geom_jitter(size=0.5, aes(color=model_clone),height=0)+ scale_color_manual(values=pal)+
  xlab('PDTs') +
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  unmute_theme+theme(legend.position="none", #axis.text.x = element_blank(), 
                     #axis.ticks.x = element_blank(),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))                   
#print(p)
  

ggsave(my_outplot, plot=pp, width=89, height=89, units="mm")

sink(my_log_f)
wilcox.test(our[our$type=="LMs", 'MR_edu'], our[our$type=="PRs", 'MR_edu'])
wilcox.test(our[our$type=="LMs", 'MR_edu'], our[our$type=="PRs", 'MR_edu'], alt='greater')
sink()