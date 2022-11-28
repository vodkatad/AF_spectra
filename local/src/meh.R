load('/scratch/trcanmed/AF_spectra/dataset/meh2.png.Rdata') # edu
our$treat <- ifelse(grepl('N', our$clone2), 'NT', 'Afatinib')
our$mt <- paste0(our$model, our$treat)
our$treat <- factor(our$treat, levels=c('NT', 'Afatinib'))
our$mt <- factor(our$mt, levels=c('CRC1430NT', 'CRC1430Afatinib','CRC1620NT', 'CRC1620Afatinib'))

#ggplot(pdata, aes(x=mt)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
ggplot(our, aes(x=mt, y=MR, color=model_clone)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
  geom_point(stat="identity", size=2)+
  geom_line(aes(group=clone))+
  scale_color_manual(values=c('#5efeef','#01e4ce','#b3fff7','#7c7c7c','#666666','#cccccc'))+ scale_x_discrete(labels=c('CRC1430NT'="NT", 'CRC1430Afatinib'="Afatinib",
                                                                               'CRC1620NT'="NT", 'CRC1620Afatinib'="Afatinib"))


load('/scratch/trcanmed/AF_spectra/dataset/meh.png.Rdata') # conte
our$treat <- ifelse(grepl('N', our$clone2), 'NT', 'Afatinib')
our$mt <- paste0(our$model, our$treat)
our$treat <- factor(our$treat, levels=c('NT', 'Afatinib'))
our$mt <- factor(our$mt, levels=c('CRC1430NT', 'CRC1430Afatinib','CRC1620NT', 'CRC1620Afatinib'))

#ggplot(pdata, aes(x=mt)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
ggplot(our, aes(x=mt, y=MR, color=model_clone)) +theme_bw()+ggtitle('MR cell counts')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
  geom_point(stat="identity", size=2)+
  geom_line(aes(group=clone))+
  scale_color_manual(values=c('#5efeef','#01e4ce','#b3fff7','#7c7c7c','#666666','#cccccc'))+ scale_x_discrete(labels=c('CRC1430NT'="NT", 'CRC1430Afatinib'="Afatinib",
                                                                                                                       'CRC1620NT'="NT", 'CRC1620Afatinib'="Afatinib"))




load('/scratch/trcanmed/AF_spectra/dataset/meh3.png.Rdata') # n muts
our$treat <- ifelse(grepl('N', our$clone2), 'NT', 'Afatinib')
our$mt <- paste0(our$model, our$treat)
our$treat <- factor(our$treat, levels=c('NT', 'Afatinib'))
our$mt <- factor(our$mt, levels=c('CRC1430NT', 'CRC1430Afatinib','CRC1620NT', 'CRC1620Afatinib'))
#ggplot(pdata, aes(x=mt)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
ggplot(our, aes(x=mt, y=MR_edu, color=model_clone)) +theme_bw()+ggtitle('Accumulated muts')+ylab('muts')+xlab('')+
  geom_point(stat="identity", size=2)+
  geom_line(aes(group=clone))+
  scale_color_manual(values=c('#5efeef','#01e4ce','#b3fff7','#7c7c7c','#666666','#cccccc'))+ scale_x_discrete(labels=c('CRC1430NT'="NT", 'CRC1430Afatinib'="Afatinib",
                                                                                                                       'CRC1620NT'="NT", 'CRC1620Afatinib'="Afatinib"))



load('/scratch/trcanmed/AF_spectra/dataset/meh4.png.Rdata') # normalized muts
our$treat <- ifelse(grepl('N', our$clone2), 'NT', 'Afatinib')
our$mt <- paste0(our$model, our$treat)
our$treat <- factor(our$treat, levels=c('NT', 'Afatinib'))
our$mt <- factor(our$mt, levels=c('CRC1430NT', 'CRC1430Afatinib','CRC1620NT', 'CRC1620Afatinib'))
our$MR <- our$MR_edu * 1000000
#ggplot(pdata, aes(x=mt)) +theme_bw()+ggtitle('MR EDU')+ylab('MR, mut/(division*bp) *10^-9')+xlab('')+
ggplot(our, aes(x=mt, y=MR, color=model_clone)) +theme_bw()+ggtitle('Accumulated muts')+ylab('muts/Mbp')+xlab('')+
  geom_point(stat="identity", size=2)+
  geom_line(aes(group=clone))+
  scale_color_manual(values=c('#5efeef','#01e4ce','#b3fff7','#7c7c7c','#666666','#cccccc'))+ scale_x_discrete(labels=c('CRC1430NT'="NT", 'CRC1430Afatinib'="Afatinib",
                                                                                                                       'CRC1620NT'="NT", 'CRC1620Afatinib'="Afatinib"))




1867  cut -f 2,7 */platypus_nobin/all.MR_baseline_ov | grep -v end > meh
1868  cut -f 2,8 */platypus_nobin/all.MR_baseline_ov | grep -v end > meh2
1870  cut -f 2,7 */platypus_nobin/all.MR_baseline_ov | grep -v end > meh
1871  cut -f 2,8 */platypus_nobin/all.MR_baseline_ov | grep -v end > meh2
1872  /scratch/trcanmed/AF_spectra/local/bin/MR_plot meh meh.png '#cc3300,#ff4000,#ff6633,#f607b9,#fb49ce,#f998e0,#9900ff,#ad33ff,#c266ff,#155d00,#239203,#2fc603,#239203,#77a003,#95c805,#95c805,#bcfc08,#bcfc08,#0829fc,#0829fc,#4a62fb,#4a62fb,#95a3fd,#95a3fd,#003399,#003399,#ff9900,#ffad33,#ff9900,#ffff00,#ffff66' /mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/theme_20.Rdata
1873  /scratch/trcanmed/AF_spectra/local/bin/MR_plot meh2 meh2.png '#cc3300,#ff4000,#ff6633,#f607b9,#fb49ce,#f998e0,#9900ff,#ad33ff,#c266ff,#155d00,#239203,#2fc603,#239203,#77a003,#95c805,#95c805,#bcfc08,#bcfc08,#0829fc,#0829fc,#4a62fb,#4a62fb,#95a3fd,#95a3fd,#003399,#003399,#ff9900,#ffad33,#ff9900,#ffff00,#ffff66' /mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/theme_20.Rdata


egrassi@godot:/scratch/trcanmed/AF_spectra/dataset$ /scratch/trcanmed/AF_spectra/local/bin/plot_EDU_conte MR_edu_disney.png '#5efeef,#2b847c,#7c7c7c,#323232' pippo pluto paperino minnie

cp NT_CRC1620/platypus_nobin/all.MR_baseline_ov paperino


egrassi@godot:/scratch/trcanmed/AF_spectra/dataset$ cat NT_CRC1430/platypus_nobin/dnds.tsv T_CRC1430/platypus_nobin/dnds.tsv NT_CRC1620/platypus_nobin/dnds.tsv  T_CRC1620/platypus_nobin/dnds.tsv | grep -v name | grep wall > paperoga_dnds.tsv
egrassi@godot:/scratch/trcanmed/AF_spectra/dataset$ vi paperoga_dnds.tsv 
egrassi@godot:/scratch/trcanmed/AF_spectra/dataset$ /scratch/trcanmed/AF_spectra/local/bin/dnds_plot_overall paperoga_dnds.tsv paperoga_dnds.png  '#5efeef,#2b847c,#7c7c7c,#323232'


#egrassi@godot:/scratch/trcanmed/AF_spectra/dataset$ cut -f 2,5 */platypus_nobin/all.MR_baseline_ov | grep -v end > meh3
#egrassi@godot:/scratch/trcanmed/AF_spectra/dataset$ cut -f 2,5,6 */platypus_nobin/all.MR_baseline_ov | grep -v end | bawk '{print $1,$2/$3}' > meh4

# mut
/scratch/trcanmed/AF_spectra/local/bin/MR_plot meh3 meh3.png '#cc3300,#ff4000,#ff6633,#f607b9,#fb49ce,#f998e0,#9900ff,#ad33ff,#c266ff,#155d00,#239203,#2fc603,#239203,#77a003,#95c805,#95c805,#bcfc08,#bcfc08,#0829fc,#0829fc,#4a62fb,#4a62fb,#95a3fd,#95a3fd,#003399,#003399,#ff9900,#ffad33,#ff9900,#ffff00,#ffff66' /mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/theme_20.Rdata


# normaliz mut
/scratch/trcanmed/AF_spectra/local/bin/MR_plot meh4 meh4.png '#cc3300,#ff4000,#ff6633,#f607b9,#fb49ce,#f998e0,#9900ff,#ad33ff,#c266ff,#155d00,#239203,#2fc603,#239203,#77a003,#95c805,#95c805,#bcfc08,#bcfc08,#0829fc,#0829fc,#4a62fb,#4a62fb,#95a3fd,#95a3fd,#003399,#003399,#ff9900,#ffad33,#ff9900,#ffff00,#ffff66' /mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/theme_20.Rdata


load('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/CRC1599PRX0A02002TUMD03000V2.fit.0.025_0.2.pdf.debug.RData')
##ff9900,#ffff00
hist(exsubcl_nohigh, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col='#ff9900', main="", ylim=c(0,30))

dev.off()
load('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/CRC1599LMX0A02001TUMD03000V2.fit.0.025_0.2.pdf.debug.RData')
##ff9900,#ffff00
hist(exsubcl_nohigh, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col='#ffff00', main="")

setwd('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/')

w <- c('CRC0053LMX0A02204TUMD07000V2.fit.0.12_0.24.pdf.debug.RData','CRC0065PRX0A01201TUMD03000V2.fit.0.12_0.24.pdf.debug.RData', 'CRC1387LMX0A02001TUMD03000V2.fit.0.12_0.24.pdf.debug.RData',
       'CRC1473PRX0B01001TUMD03000V2.fit.0.12_0.24.pdf.debug.RData', 'CRC0053PRX0A01201TUMD05000V2.fit.0.12_0.24.pdf.debug.RData','CRC1331LMX0A02005TUMD03000V2.fit.0.12_0.24.pdf.debug.RData',
       'CRC1387PRX0B01001TUMD05000V2.fit.0.12_0.24.pdf.debug.RData','CRC1599LMX0A02001TUMD03000V2.fit.0.12_0.24.pdf.debug.RData', 'CRC0065LMX0B02205TUMD02000V2.fit.0.12_0.24.pdf.debug.RData', 
       'CRC1331PRX0A01001TUMD05000V2.fit.0.12_0.24.pdf.debug.RData','CRC1473LMX0B02002TUMD03000V2.fit.0.12_0.24.pdf.debug.RData')
w <- w[order(w)]
i <- 0
xlim <- c(30, 30, 40, 40, 50, 50, 50, 50, 40, 40)
j <- 1
for (rd in w) {
  if (i == 0) {
    col = '#595959'
    i = 1
  } else {
    col = '#adacac'
    i = 0
  }
  load(rd)
  model <- substr(rd, 0, 10)
  hist(exsubcl_nohigh, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col=col, ylim=c(0,lim[j]), main=model)
  j <- j + 1
}

