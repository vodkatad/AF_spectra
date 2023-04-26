library(ggplot2)
library(showtext)
size <- as.numeric(snakemake@wildcards[['size']])

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

# Da Marti e https://www.christophenicault.com/post/understand_size_dimension_ggplot2/
showtext_opts(dpi = 300) 
# since we are not changing fonts in the end cause myriad end up not being text object I'm not sure it's needed
showtext_auto(enable = TRUE)

textSize <- size
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
	panel.background = element_blank()
)


mute_theme <- theme_bw() +
theme(
	text = element_blank()
)

save.image(snakemake@output[['Rimage']])

#pdata <- data.frame(x=rnorm(100), y=rnorm(100))
#ggplot(data=pdata, aes(x, y))+ggtitle('example')+geom_point()+unmute_theme
#ggsave('test.png')
# width = 8, height = 7, units = "in", dpi = 300

#egrassi@qui:~/Dropbox/work/biobanca/figures/3b/def$ perl -pne "s/textLength='.+px'//; " < pdo_vs_msk.svg > pdo_vs_msk_notexlen.svg
