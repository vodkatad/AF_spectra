data_f  <- snakemake@input[['data']]
colors <- snakemake@input[['palette']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]

theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(ggplot2)
load(theme)

palette_df <- readRDS(colors)
pal <- palette_df$palette
names(pal) <- palette_df$model

load(data_f)
our$treat <- ifelse(grepl('N', our$clone2), 'NT', 'Afatinib')
our$mt <- paste0(our$model, our$treat)
our$treat <- factor(our$treat, levels=c('NT', 'Afatinib'))
our$mt <- factor(our$mt, levels=c('CRC1430NT', 'CRC1430Afatinib','CRC1620NT', 'CRC1620Afatinib'))

y_breaks <- guess_ticks(our$MR)


p <- ggplot(our, aes(x=mt, y=MR, color=model_clone)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
  geom_point(stat="identity", size=2)+
  geom_line(aes(group=clone))+
  scale_color_manual(values=pal)+
  scale_y_continuous(breaks=y_breaks,limits=c(0, max(y_breaks)),expand = c(0, 0))+# + ylim(min(y_breaks),max(y_breaks))+
  scale_x_discrete(labels=c('CRC1430NT'="NT", 'CRC1430Afatinib'="Afatinib", 'CRC1620NT'="NT", 'CRC1620Afatinib'="Afatinib"))+
  unmute_theme+theme(legend.position="none", #axis.text.x=element_text(angle = -90, hjust = 0),
                     legend.spacing.y = unit(0.15, "mm")) + guides(col=guide_legend(nrow=length(pal), keyheight=unit(0.01, "mm")))   


ggsave(outplot, plot=p, width=89, height=89, units="mm")

save.image(paste0(outplot, '.Rdata'))
