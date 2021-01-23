setwd('/scratch/trcanmed/AF_spectra/dataset/')
load('pluto.Rdata')

wp <- function(pair,data) {
  pair <- unlist(pair)
  p1 <- data[data$model==pair[1],'MR']
  p2 <- data[data$model==pair[2],'MR']
  wilcox.test(p1, p2)$p.value
}

pairs <- expand.grid(unique(our$model),unique(our$model))

wilcoxpval <- matrix(apply(pairs, 1, wp, our),nrow=length(unique(our$model)))
colnames(wilcoxpval) <- unique(our$model)
rownames(wilcoxpval) <- unique(our$model)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
"cyan", "#007FFF", "blue", "#00007F"))
corrplot(-log10(wilcoxpval), method = "circle", is.corr=FALSE, col=col1(100))

lpval <- -log10(wilcoxpval)
minv <- min(lpval)
maxv <- max(lpval)
neutral_value <- 2
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.1),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.1))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
                "snow1", "snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

#corrplot(lpval, method = "circle", is.corr=FALSE, col=my_palette)

corrplot.mixed(lpval, is.corr=FALSE, lower.col = "black", number.cex = .7, col.upper = my_palette)

pheatmap(lpval, color=my_palette, cluster_cols=F, cluster_rows=F)

lpval <- -log10(wilcoxpval)
lpval2 <- lpval
lpval2[lower.tri(lpval)] <- round(lpval2[lower.tri(lpval)], 2)
lpval2[upper.tri(lpval)] <- ''
lpval[lower.tri(lpval)] <-   2.088389
pheatmap(lpval, color=my_palette, cluster_cols=F, cluster_rows=F, display_numbers=lpval2)

diffs <- function(pair,data) {
  pair <- unlist(pair)
  p1 <- data[rownames(data)=='mean', colnames(data)==pair[1]]
  p2 <- data[rownames(data)=='mean', colnames(data)==pair[2]]
  abs(p1-p2)
}

meandiff <- matrix(apply(pairs, 1, diffs, ic_clones),nrow=length(unique(our$model)))
colnames(meandiff) <- unique(our$model)
rownames(meandiff) <- unique(our$model)

lpval <- meandiff
minv <- min(lpval)
maxv <- max(lpval)

bk <- c(seq(minv-0.1,maxv+0.1,by=0.001))
my_palette <- c(colorRampPalette(colors = c("snow1","darkblue"))(n = length(bk)),'snow1')


lpval2 <- lpval
lpval2[lower.tri(lpval)] <- round(lpval2[lower.tri(lpval)], 2)
lpval2[upper.tri(lpval)] <- ''
lpval[lower.tri(lpval)] <-   10
pheatmap(lpval, color=my_palette, cluster_cols=F, cluster_rows=F, display_numbers=lpval2)
pheatmap(lpval, color=rev(brewer.pal(length(my_palette), 'Greens')), cluster_cols=F, cluster_rows=F, display_numbers=lpval2)


###
lmp <- function (sm) {
  #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  #f <- summary(modelobject)$fstatistic
  f <- sm$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


lm_model <- function(model,data) {
  da <- data[data$model==model,]
  sm <- summary(lm(data=da, formula(MR~clone)))
  lmp(sm)
}

lm_s_model <- function(model, data) {
  da <- data[data$model==model,]
  da$MR <- sample(da$MR)
  sm <- summary(lm(data=da, formula(MR~clone)))
  lmp(sm)
}

###
vivo <- read.table('/scratch/trcanmed/AF_spectra/dataset/vivo_gained_SNV', header=F, sep="\t", stringsAsFactors=FALSE)
vitro <- read.table('/scratch/trcanmed/AF_spectra/dataset/MR_edu_SNV', header=F, sep="\t", stringsAsFactors=FALSE)
colnames(vivo) <- c('sample','n')
colnames(vitro) <- c('sample','MR')

add_info <- function(our) {
  our$model <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
  our$clone <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][2]; return(y[1])})
  our$clone2 <- sapply(our$sample, function(x) {y<-strsplit(x, '-')[[1]][4]; return(y[1])})
  our$modelclone <- paste0(our$model, '-', our$clone)
  return(our)
}

vivo <- add_info(vivo)
vitro <- add_info(vitro)

# now need to average on modelclone


vitro_mean <- as.data.frame(sapply(unique(vitro$modelclone), function(x) { mean(vitro[vitro$modelclone==x,'MR']) }))
vivo_mean <- as.data.frame(sapply(unique(vivo$modelclone), function(x) { mean(vivo[vivo$modelclone==x,'n']) }))
colnames(vivo_mean)[1] <- 'n'
colnames(vitro_mean)[1] <- 'MR'
m <- merge(vitro_mean, vivo_mean, by="row.names")
plot(m$MR, m$n)
cor.test(m$MR, m$n)
cor.test(m$MR, m$n, method="spearman")
m$model <- sapply(m$Row.names, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})
ggplot(m, aes(x=MR, y=n, color=model)) +  geom_point()


##

vitro_mean <- as.data.frame(sapply(unique(vitro$model), function(x) { mean(vitro[vitro$model==x,'MR']) }))
vivo_mean <- as.data.frame(sapply(unique(vivo$model), function(x) { mean(vivo[vivo$model==x,'n']) }))
colnames(vivo_mean)[1] <- 'n'
colnames(vitro_mean)[1] <- 'MR'
m <- merge(vitro_mean, vivo_mean, by="row.names")
plot(m$MR, m$n)
cor.test(m$MR, m$n)
cor.test(m$MR, m$n, method="spearman")
m$model <- sapply(m$Row.names, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})

ggplot(m, aes(x=MR, y=n, color=model)) +  geom_point(size=3)+theme_bw()+scale_color_manual(values=c(TODO))



