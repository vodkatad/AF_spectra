load('/scratch/trcanmed/AF_spectra/dataset/meh.png.Rdata')
our$treat <- ifelse(grepl('N', our$clone2), 'NT', 'Afatinib')
our$mt <- paste0(our$model, our$treat)
#ggplot(pdata, aes(x=mt)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
geom_point(data=our, aes(x=mt, y=MR, color=treat), stat="identity", size=2, position=position_dodge(0.2))
ggplot(pdata, aes(x=mt)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
geom_point(data=our, aes(x=mt, y=MR, color=treat), stat="identity", size=2, position=position_dodge(0.2))
ggplot(pdata, aes(x=mt)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
geom_point(data=our, aes(x=mt, y=MR, color=treat), stat="identity", size=2
)
