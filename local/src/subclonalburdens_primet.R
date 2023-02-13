library(ggplot2)
data <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/mutect_paired/subcl_mut_burden.tsv', sep= "\t", header=T)

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
ggsave('/home/egrassi/primet_subcl.pdf')


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

ggplot(data=d3, aes(x=mp, y=burden, fill=mp))+
  geom_violin()+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Original Tumor')+ylab('Subclonal mutational burden (/ Mbps)')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title="")) + stat_summary(fun.data=data_summary)
ggsave('/home/egrassi/primet_subcl_vio.pdf')

met <- data[data$mp == "LMX",]
pri <- data[data$mp == "PRX",]
if (!all(met$smodel==pri$smodel)) {
  stop('llama! Qualquadra non cosa  pri-met pairs')
}
t.test(met$burden, pri$burden, alternative="greater", paired=TRUE)
wilcox.test(met$burden, pri$burden, alternative="greater", paired=TRUE)




ggplot(data=d3, aes(x=mp, y=burden, fill=mp))+
  geom_violin()+theme_bw()+
  geom_point()+geom_line(aes(group=smodel, color=as.factor(smodel)))+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Original Tumor')+ylab('Subclonal mutational burden (/ Mbps)')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))
ggsave('/home/egrassi/primet_subcl_bp.pdf')

### slope

data <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/bestbet_0.12_0.24.tsv', sep="\t", header=TRUE)

data$lmodel <- substr(rownames(data), 0, 10)
data$smodel <- substr(rownames(data), 0, 7)
data$mp <- substr(rownames(data), 8, 10)
data$mp <- factor(data$mp, levels=c('PRX', 'LMX'))

fit_r2 <- data
fit_keep <- fit_r2[fit_r2$r > 0.90 & fit_r2$subcl > 10,]

pairs <- as.data.frame(table(fit_keep$smodel))
with_pair <- fit_keep[fit_keep$smodel %in% pairs[pairs$Freq == 2,'Var1'],]
fit_r2 <- with_pair
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}


ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

LEVEL <- 0.99
ic_clones <- sapply(unique(fit_r2$mp), function(x) { confidence_interval(fit_r2[fit_r2$mp==x,'intercept'], LEVEL) })
colnames(ic_clones) <- unique(fit_r2$mp)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

pdata$model <- factor(pdata$model, levels=c('PRX', 'LMX'))
fit_r2$model <- factor(fit_r2$model, levels=c('PRX', 'LMX'))
ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=5) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+ctheme+
  geom_jitter(data=fit_r2, aes(x=mp, y=intercept, color=mp), size=4, shape=18, width=0.1)+
  scale_color_manual(values=c('#adacac', '#595959'))+
  ggtitle('Estimated MR')+ylab('μ/β')+xlab('')+
  geom_signif(data=fit_r2, mapping=aes(x=mp, y=intercept), 
              comparisons = list(c("LMX", "PRX")), test="t.test", test.args=list(alternative = "greater", paired=TRUE))

pd <- position_dodge(width=0.2)
  
#https://stackoverflow.com/questions/39533456/how-to-jitter-both-geom-line-and-geom-point-by-the-same-magnitude
ggplot(data=fit_r2, aes(x=mp, y=intercept, color=mp, group=smodel)) +
  geom_jitter(data=fit_r2, aes(x=mp, y=intercept, color=mp), size=4, shape=18, position=pd)+
  geom_line(data=fit_r2, aes(group=smodel), position=pd, color="lightgrey", linetype = "dashed")+
  geom_point(data=pdata, aes(x=model, y=mean, group=NULL, color=NULL), stat="identity", shape=1, size=5) +
  geom_segment(data=pdata, aes(y=lower, yend=upper, x=model, xend=model, group=NULL, color=NULL), size=0.6) + ctheme+
  scale_color_manual(values=c('#adacac', '#595959'))+
  ggtitle('Estimated MR')+ylab('μ/β')+xlab('')+
  geom_signif(data=fit_r2, mapping=aes(x=mp, y=intercept), 
              comparisons = list(c("LMX", "PRX")), test="t.test", test.args=list(alternative = "greater", paired=TRUE))



met <- fit_r2[fit_r2$mp == "LMX",]
pri <- fit_r2[fit_r2$mp == "PRX",]
if (!all(met$smodel==pri$smodel)) {
  stop('llama! Qualquadra non cosa in pri-met pairs')
}
ti <- t.test(met$intercept, pri$intercept, alternative="greater", paired=TRUE)

ggsave('slopes_filtered_lines.svg', height=5.25, width=5.25, units="in")

#######33 signatures tries
data <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/signatures/all_cosmic_fit.tsv', sep="\t", header=T, row.names=1)
sdata <- data[-nrow(data),]

pdata<-data.frame(row.names=colnames(sdata), ratio=unlist((sdata[8,]+0.00001)/(sdata[1,]+0.00001)))
pdata$lmodel <- substr(rownames(pdata), 0, 10)
pdata$smodel <- substr(rownames(pdata), 0, 7)
pdata$mp <- substr(rownames(pdata), 8, 10)
pdata$mp <- factor(pdata$mp, levels=c('PRX', 'LMX'))

LEVEL <- 0.99
ic_clones <- sapply(unique(pdata$mp), function(x) { confidence_interval(pdata[pdata$mp==x,'ratio'], LEVEL) })
colnames(ic_clones) <- unique(pdata$mp)
sumdata <- as.data.frame(t(ic_clones))
sumdata$model <- rownames(sumdata)
sumdata$model <- factor(sumdata$model, levels=c('PRX', 'LMX'))


ggplot(sumdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=5) +ctheme+
  geom_jitter(data=pdata, aes(x=mp, y=ratio, color=mp), size=4, shape=18, width=0.1)+
  scale_color_manual(values=c('#adacac', '#595959'))+
  ggtitle('Signature 8/1')+ylab('Exposures ratio')+xlab('')+
  geom_signif(data=pdata, mapping=aes(x=mp, y=ratio), 
              comparisons = list(c("LMX", "PRX")), test="t.test", test.args=list(alternative = "greater", paired=TRUE))


pdata<-data.frame(row.names=colnames(sdata), S1=unlist(sdata[1,]), S8=unlist(sdata[8,]))
pdata$lmodel <- substr(rownames(pdata), 0, 10)
pdata$smodel <- substr(rownames(pdata), 0, 7)
pdata$mp <- substr(rownames(pdata), 8, 10)
pdata$mp <- factor(pdata$mp, levels=c('PRX', 'LMX'))

library(reshape)
me <- melt(pdata, id="lmodel", values=c('S1', 'S8'))

me$mp <- substr(rownames(me), 8, 10)
me$mp <- factor(me$mp, levels=c('PRX', 'LMX'))

ggplot(data=pdata, aes(x=smodel, y=S1, fill=mp)) + theme_bw()+
  geom_col(position='dodge')+
  scale_fill_manual(values=c('#adacac', '#595959'))
  

ggplot(data=pdata, aes(x=smodel, y=S8, fill=mp)) + theme_bw()+
  geom_col(position='dodge')+
  scale_fill_manual(values=c('#adacac', '#595959'))


###### conta subclonale pura

data <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/mutect_paired/subcl_mut_burden.tsv', sep= "\t", header=T)

rownames(data) <- data$X
data$lmodel <- substr(rownames(data), 0, 10)
data$smodel <- substr(rownames(data), 0, 7)
data$mp <- substr(rownames(data), 8, 10)
data$mp <- factor(data$mp, levels=c('PRX', 'LMX'))


fit_r2 <- data
#fit_keep <- fit_r2[fit_r2$r > 0.90 & fit_r2$subcl > 10,]

#pairs <- as.data.frame(table(fit_keep$smodel))
#with_pair <- fit_keep[fit_keep$smodel %in% pairs[pairs$Freq == 2,'Var1'],]

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, "mean" = vec_mean)
  return(result)
}


ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                           axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                           plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.position='none'
)

LEVEL <- 0.99
ic_clones <- sapply(unique(fit_r2$mp), function(x) { confidence_interval(fit_r2[fit_r2$mp==x,'totmut'], LEVEL) })
colnames(ic_clones) <- unique(fit_r2$mp)
pdata <- as.data.frame(t(ic_clones))
pdata$model <- rownames(pdata)

pdata$model <- factor(pdata$model, levels=c('PRX', 'LMX'))
fit_r2$mp <- factor(fit_r2$mp, levels=c('PRX', 'LMX'))
ggplot(pdata, aes(x=model, y=mean)) +  geom_point(stat="identity", shape=1, size=5) +
  geom_segment(aes(y=lower, yend=upper, x=model, xend=model), size=0.6)+ctheme+
  geom_jitter(data=fit_r2, aes(x=mp, y=totmut, color=mp), size=4, shape=18, width=0.1)+
  scale_color_manual(values=c('#adacac', '#595959'))+
  ggtitle('Subclonal mutations')+ylab('N.')+xlab('')+
  geom_signif(data=fit_r2, mapping=aes(x=mp, y=totmut), 
              comparisons = list(c("LMX", "PRX")), test="t.test", test.args=list(alternative = "greater", paired=TRUE))

met <- fit_r2[fit_r2$mp == "LMX",]
pri <- fit_r2[fit_r2$mp == "PRX",]
if (!all(met$smodel==pri$smodel)) {
  stop('llama! Qualquadra non cosa in pri-met pairs')
}
ti <- t.test(met$totmut, pri$totmut, alternative="greater", paired=TRUE)

ggsave('subclonal_n.pdf', height=5.25, width=5.25, units="in")
  