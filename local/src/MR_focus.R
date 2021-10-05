
load('/scratch/trcanmed/AF_spectra/dataset/MR_edu_SNV.png.Rdata')

w <- c('CRC1078','CRC1307','CRC1599LM','CRC1599PR')
pdata <- pdata[pdata$model %in% w,]
grepps <- unlist(sapply(w, function(x) { our$model[grepl(x, our$model,)] }))
our <- our[our$model %in% grepps,]

cbPalette <- c('#155d00','#239203','#2fc603','#77a003','#95c805','#95c805','#bcfc08','#bcfc08','#ffff00','#ffff66','#ff9900','#ffad33')

pdata$model <- factor(pdata$model, levels=c('CRC1078','CRC1307','CRC1599PR','CRC1599LM'))

library(ggsignif)

ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
  geom_point(data=our, aes(x=model, y=MR, color=model_clone, shape=time), stat="identity", size=4, position=position_dodge(0.2))+
  ctheme+scale_color_manual(values=cbPalette)+scale_shape_manual(values=c(18,20))+geom_signif(data=our, aes(x=model, y=MR),
    comparisons = list(c("CRC1078", "CRC1307"), c('CRC1599PR','CRC1599LM')),
    map_signif_level = TRUE, manual=TRUE
    #test="t.test"
  )

wilcox.test(our[grepl('CRC1078',our$sample),'MR_edu'], our[grepl('CRC1307',our$sample),'MR_edu'])

wilcox.test(our[grepl('CRC1599LM',our$sample),'MR_edu'], our[grepl('CRC1599PR',our$sample),'MR_edu'])


load('/scratch/trcanmed/AF_spectra/dataset/MR_edu_SNV.png.Rdata')

ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15), axis.title.x = element_text(size=20),
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

w <- c('CRC1502')
pdata <- pdata[pdata$model %in% w,]
grepps <- unlist(sapply(w, function(x) { our$model[grepl(x, our$model,)] }))
our <- our[our$model %in% grepps,]
cbPalette <- c('#0829fc','#0829fc','#4a62fb','#4a62fb','#95a3fd','#95a3fd','#003399','#003399')

ggplot(our, aes(x=time, y=MR, shape=time, color=model_clone)) +  geom_point(stat="identity", size=3, position=position_dodge(0.2)) +
  ctheme+scale_shape_manual(values=c(18,20))+scale_color_manual(values=cbPalette)+geom_signif(color='black',comparisons = list(c("1", "2")),
                                                                                              map_signif_level = TRUE, 
                                                                                              #test="t.test"
  )

wilcox.test(our[grepl('1',our$time),'MR_edu'], our[grepl('2',our$time),'MR_edu'])
