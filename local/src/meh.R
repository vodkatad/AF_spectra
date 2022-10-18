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