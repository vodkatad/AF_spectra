

library(vcfR)
library(ggplot2)
basedir <- '/scratch/trcanmed/AF_spectra/local/share/data/cetuxi_treat_paneth/'
CRC0542 <- read.vcfR(paste0(basedir, '/anot_CRC542_all_platypus_filtered.vcf'))
CRC0252 <- read.vcfR(paste0(basedir,'/anot_CRC252_all_platypus_filtered.vcf'))
map <- read.table(paste0(basedir, '/samples.txt'), sep="\t", header=TRUE)



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
## N of clonal mutations


#keep only PASS?
#no in the suppl they say that they keep everything in ther:
# data@rotpunkt:~/work/stash/tesi_ire/sottoriva_scatter$ grep -v "^#" anot_CRC542_all_platypus_filtered.vcf| cut -f7  |sort | uniq
# alleleBias
# HapScore
# PASS
# Q20
# Q20;alleleBias
# Q20;QD
# Q20;QD;alleleBias
# QD
# QD;alleleBias
# SC;alleleBias
# SC;QD;alleleBias
# data@rotpunkt:~/work/stash/tesi_ire/sottoriva_scatter$ grep -v "^#" anot_CRC252_all_platypus_filtered.vcf| cut -f7  |sort | uniq
# alleleBias
# HapScore
# HapScore;alleleBias
# PASS
# Q20
# Q20;alleleBias
# Q20;QD
# Q20;QD;alleleBias
# QD
# QD;alleleBias
# SC;alleleBias


get_vaf <- function(vcf) {
  depth <- extract.gt(vcf, element = "NR", as.numeric = TRUE)
  # How many reads support the variant?
  variant <- extract.gt(vcf, element = "NV", as.numeric = TRUE)
  # Calculate the vaf
  vafs <- variant / depth
  fix <- as.data.frame(getFIX(vcf))
  res <- data.frame(chr=fix$CHROM, pos=fix$POS, ref=fix$REF, alt=fix$ALT, vaf=vafs)
}

CRC0542_vaf <- as.data.frame(get_vaf(CRC0542))
CRC0252_vaf <- as.data.frame(get_vaf(CRC0252))

CRC0252_vaf$normal <- NULL
CRC0542_vaf$normal <- NULL
colnames(CRC0252_vaf) <- gsub('_merged_duplicates','', colnames(CRC0252_vaf))
colnames(CRC0542_vaf) <- gsub('_merged_duplicates','', colnames(CRC0542_vaf))

colnames(CRC0252_vaf) <- gsub('vaf.','', colnames(CRC0252_vaf), fixed=TRUE)
colnames(CRC0542_vaf) <- gsub('vaf.','', colnames(CRC0542_vaf), fixed=TRUE)


clonal_strict <- function(vaf) {
  vaf[apply(vaf, 1, function(x) {all(x > 0.3) }),]
}

clonal_lax <- function(vaf) {
  vaf[apply(vaf, 1, function(x) {any(x > 0.3) }),]
}


fixnames <- function(vaf, map, name) {
  #map <- map[,c('hidden','description')]
  map$description <- as.character(map$description)
  map$description <- gsub('WES ', '', map$description)
  map$description <- gsub(paste0(' ', name), '', map$description)
  map$hidden <- as.character(map$hidden)
  vaf_keep <- vaf[,colnames(vaf) %in% map$hidden]
  vaf_add <- vaf[,!colnames(vaf) %in% map$hidden]
  map <- map[match(colnames(vaf_keep), map$hidden),]
  colnames(vaf_keep) <- map$description
  return(cbind(vaf_keep, vaf_add))
}


mat252 <- fixnames(CRC0252_vaf, map, 'CRC0542') #error in the samples.gnumeric file in description TODO RICAPIRE
mat542 <- fixnames(CRC0542_vaf, map, 'CRC0542')

summary(mat252$normal) # all 0
summary(mat542$normal) # all 0

cols <- c('chr', 'pos', 'ref', 'alt')
id252 <- mat252[, cols]
id542 <- mat542[, cols]

cols <- c(cols, 'normal')
mat252[, cols] <- NULL
mat542[, cols] <- NULL

n_clonal252 <- apply(mat252, 2, function(x) { sum(x>=0.3) })
n_clonal542 <- apply(mat542, 2, function(x) { sum(x>=0.3) })


plotclon <- function(n_clonal,name, ylab) {
  data <- data.frame(model=rep(name,length(n_clonal)), n_clonal = n_clonal, treat=names(n_clonal))
  data$treat <- as.factor(data$treat)
  data$treat <- factor(data$treat, levels=c('Placebo','Cetux','Release'))
  data <- data[order(data$treat),]
  data$x <- seq(1,nrow(data))
  print(ggplot(data=data, aes(x=x, y=n_clonal, fill=treat))+geom_col()+current_theme+ggtitle(name)+xlab('')+ylab(ylab)+ylim(0, 350))
  #ggsave(paste0('n_clonal', name, '.svg'))
  #write.xlsx(data, paste0('n_clonal', name, '.xlsx'))
}

plotclon(n_clonal252,'CRC0252', '# clonal muts')
plotclon(n_clonal542,'CRC0542', '# clonal muts')

n_sclonal252 <- apply(mat252, 2, function(x) { sum(x<=0.3 & x > 0.1) })
n_sclonal542 <- apply(mat542, 2, function(x) { sum(x<=0.3 & x > 0.1) })

plotclon(n_sclonal252,'CRC0252', '# subclonal muts')
plotclon(n_sclonal542,'CRC0542', '# subclonal muts')

setwd('/scratch/trcanmed/tmp/au')

for (i in seq(1, ncol(mat252))) {
  name <- paste0('CRC0252_', colnames(mat252)[i], '_', i)
  data <- cbind(name, id252, mat252[,i])
  data <- data[data[,6] > 0.1,]
  write.table(data, file=paste0('dndsin_', name, '.tsv'), sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
}
#egrassi@godot:/scratch/trcanmed/tmp/au$ /scratch/trcanmed/AF_spectra/local/bin/dnds dndsin_CRC0252_Cetux.tsv o i /scratch/trcanmed/AF_spectra/local/share/data/RefCDS_human_GRCh38.p12.rda 
#for i in dndsin*; do /scratch/trcanmed/AF_spectra/local/bin/dnds $i dndsout_$i tmp /scratch/trcanmed/AF_spectra/local/share/data/RefCDS_human_GRCh38.p12.rda; done;



d <- read.table('dndsi.tsv', sep="\t", header=F, stringsAsFactors = F)
d[,2] <- NULL
d[,2] <- NULL
colnames(d) <- c('model','estimate', 'upper', 'lower')
d$model <- gsub('dndsout_dndsin_', '', d$model)
d$model <- gsub('.tsv', '', d$model, fixed=T)



ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)


ggplot(d, aes(x=model, y=estimate, color=model)) +  geom_point(stat="identity", size=2) +
  geom_errorbar(aes(ymin=lower, ymax=upper, x=model, width=0.1, color=model), size=0.2)+theme_bw()+ggtitle('dN/dS nonsyn+truncating')+ylab('ML estimate')+xlab('')+
  ctheme
#### scatters
average <- function(data) {
  data_ce <- data[, grepl('Cetux', colnames(data))]
  data_ve <- data[, grepl('Placebo', colnames(data))]
  data_re <- data[, grepl('Release', colnames(data))]
  me_ce <- as.data.frame(rowMeans(data_ce))
  me_ve <- as.data.frame(rowMeans(data_ve))
  me_re <- as.data.frame(rowMeans(data_re))
  colnames(me_ce) <- ('Cetux_VAF')
  colnames(me_ve) <- ('Placebo_VAF')
  colnames(me_re) <- ('Release_VAF')
  m <- merge(me_ce, me_ve, by="row.names")
  rownames(m) <- m$Row.names
  m$Row.names <- NULL
  m2 <- merge(m, me_re, by="row.names")
  return(m2)
}

nave252 <- average(mat252)
nave542 <- average(mat542)


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#https://slowkow.com/notes/ggplot2-color-by-density/
compare <- function(x, y, log, nx, ny, title) {
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y, density=get_density(x,y, n=100))
  ggplot(d, aes(x=x, y=y, color=density)) +geom_point()+current_theme+xlab(nx)+ylab(ny)+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+scale_color_viridis_c()+ggtitle(title)
  #ggsave(paste0(title, "_",nx,"_",ny,'.svg'), width=10, height=10, units="in")
}

compare(nave542$Placebo_VAF, nave542$Cetux_VAF, FALSE,'Placebo - mean VAF','Cetux - mean VAF', 'CRC0542')
compare(nave542$Placebo_VAF, nave542$Release_VAF, FALSE,'Placebo - mean VAF','Release - mean VAF', 'CRC0542')

compare(nave252$Placebo_VAF, nave252$Cetux_VAF, FALSE,'Placebo - mean VAF','Cetux - mean VAF', 'CRC0252')
compare(nave252$Placebo_VAF, nave252$Release_VAF, FALSE,'Placebo - mean VAF','Release - mean VAF', 'CRC0252')

#write.xlsx(nave542, 'scatter_WES_542.xlsx')
#write.xlsx(nave252, 'scatter_WES_252.xlsx')
