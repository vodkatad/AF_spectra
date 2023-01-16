library(ggplot2)
data <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/mutect_paired/clonal_mut_burden.tsv', sep= "\t", header=T)

rownames(data) <- data$X
data$lmodel <- substr(rownames(data), 0, 10)
data$smodel <- substr(rownames(data), 0, 7)
data$mp <- substr(rownames(data), 8, 10)
data$mp <- factor(data$mp, levels=c('PRX', 'LMX'))

lar <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'burden'] - sub[sub$mp=="PRX",'burden']
}


meslo <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'burden']
}


dd2 <- as.data.frame(sapply(unique(data$smodel), lar, data))
dd2 <- dd2[order(dd2[,1]),, drop=FALSE]

dd3 <- as.data.frame(sapply(unique(data$smodel), meslo, data))
dd3 <- dd3[order(dd3[,1]),, drop=FALSE]



d3 <- merge(dd2, data, by.x="row.names", by.y="smodel")
colnames(d3)[2] <- 'o'
colnames(d3)[1] <- 'smodel'
ggplot(data=d3, aes(x=reorder(smodel,o), y=burden, fill=mp))+
  geom_col(position="dodge")+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Patient')+ylab('Subclonal mutational burden (/ Mbps)')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))


met <- data[data$mp == "LMX",]
pri <- data[data$mp == "PRX",]
if (!all(met$smodel==pri$smodel)) {
  stop('llama! Qualquadra non cosa  pri-met pairs')
}
t.test(met$burden, pri$burden, alternative="greater", paired=TRUE)
wilcox.test(met$burden, pri$burden, alternative="greater", paired=TRUE)
ggsave('/home/egrassi/primet_subcl.pdf')

### slopes

d <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/bestbet_0.1_0.2.tsv', sep="\t", header=TRUE)

d$lmodel <- substr(rownames(d), 0, 10)


d$smodel <- substr(rownames(d), 0, 7)
d$mp <- substr(rownames(d), 8, 10)

d$mp <- factor(d$mp, levels=c('PRX', 'LMX'))

d2 <- d[d$r > 0.90 & d$subcl > 10,]

dd <- as.data.frame(table(d2$smodel))
d2 <- d2[d2$smodel %in% dd[dd$Freq == 2,'Var1'],]

ggplot(data=d2, aes(x=smodel, y=intercept, fill=mp))+
  geom_col(position="dodge")+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Patient')+ylab('Slope μ/β')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))

lar <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'intercept'] - sub[sub$mp=="PRX",'intercept']
}


meslo <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'intercept']
}


dd2 <- as.data.frame(sapply(unique(d2$smodel), lar, d2))
dd2 <- dd2[order(dd2[,1]),, drop=FALSE]

dd3 <- as.data.frame(sapply(unique(d2$smodel), meslo, d2))
dd3 <- dd3[order(dd3[,1]),, drop=FALSE]



d3 <- merge(dd2, d2, by.x="row.names", by.y="smodel")
colnames(d3)[2] <- 'o'
colnames(d3)[1] <- 'smodel'
ggplot(data=d3, aes(x=reorder(smodel,o), y=intercept, fill=mp))+
  geom_col(position="dodge")+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Patient')+ylab('Slope μ/β')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))