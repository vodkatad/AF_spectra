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

