library(ggdendro)

data <- read.table('/scratch/trcanmed/AF_spectra/dataset/pc_naive_trees/binary_matrix.tsv.gz', sep="\t", header=TRUE, quote="")
rownames(data) <- paste0(data[,'CHROM'],":", data[,'POS'], ":", data[,'REF'], ":", data[,'ALT'])
data[,'CHROM'] <- NULL
data[,'POS'] <- NULL
data[,'REF'] <- NULL
data[,'ALT'] <- NULL

d <- data[, !grepl('NMH', colnames(data))] 
d <- d[, !grepl('NLH', colnames(d))] 

not_private <- rowSums(d) > 1
dall <- d[not_private,] 
#> dim(d); dim(dall)
#[1] 445647    113
#[1] 126062    113
model <- hclust(dist(t(dall)), "ward.D2")

#dendro <- as.dendrogram(model)
#ddata <- dendro_data(dendro, type="rectangle")
#ggdendrogram(ddata)

#https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
dendro_data_k <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}

plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {
  
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  
  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  
  p
}

hcdata <- dendro_data_k(model, 5)

p <- plot_ggdendro(hcdata,
                   direction   = "lr",
                   expand.y    = 0.2)


cbPalette2 <- c('black',"#ff5733", 
                #"#9d01fc"
                "#f607b9",
                "#155d00",
                "#77a003",
                "#0829fc"
)
p+scale_color_manual(values=cbPalette2)+theme_bw()


###
pca <- prcomp(t(dall))

percentVar <- pca$sdev^2/sum(pca$sdev^2)
pcs <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample = rownames(pca$x))
pcs$model <- (unlist(lapply(strsplit(as.character(pcs$sample),'.', fixed=TRUE), function(x){ x[1] })))
pcs$model2 <- substr(pcs$model, 0,7)
ggplot(data = pcs, aes_string(x = "PC1", y = "PC2", color = "model2")) +geom_point(size = 1) + xlab(paste0("PC1", ': ', round(percentVar[1] *100), "% variance")) + ylab(paste0("PC2", ":", round(percentVar[2] *100), "% variance")) + coord_fixed()+scale_color_manual(values=cbPalette2[-1])+theme_bw()

# no msi

dallnomsi <- dall[, !grepl('CRC0282',colnames(dall) )]
not_private <- rowSums(dallnomsi) > 1
dallnomsi <- dallnomsi[not_private,] 

pca <- prcomp(t(dallnomsi))

percentVar <- pca$sdev^2/sum(pca$sdev^2)
pcs <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample = rownames(pca$x))
pcs$model <- (unlist(lapply(strsplit(as.character(pcs$sample),'.', fixed=TRUE), function(x){ x[1] })))
pcs$model2 <- substr(pcs$model, 0,7)

cbPalette3 <- c(
                "#f607b9",
                "#155d00",
                "#77a003",
                "#0829fc"
)
ggplot(data = pcs, aes_string(x = "PC1", y = "PC2", color = "model2")) +geom_point(size = 1) + xlab(paste0("PC1", ': ', round(percentVar[1] *100), "% variance")) + ylab(paste0("PC2", ":", round(percentVar[2] *100), "% variance")) + coord_fixed()+scale_color_manual(values=cbPalette3)+theme_bw()

###

cbPalette4 <- c("#ff5733", "#ff7433", "#f607b9","#fb49ce","#155d00","#239203","#77a003","#95c805","#0829fc","#4a62fb")
  
ggplot(data = pcs, aes_string(x = "PC1", y = "PC2", color = "model")) +geom_point(size = 1) + xlab(paste0("PC1", ': ', round(percentVar[1] *100), "% variance")) + ylab(paste0("PC2", ":", round(percentVar[2] *100), "% variance")) + coord_fixed()+theme_bw()+scale_color_manual(values=cbPalette4)

####
pcs$time <- unlist(lapply(strsplit(as.character(pcs$sample),'.', fixed=T), function(x){ if (length(x)>2) {x[3]} else {'B'} }))
pcs$env <- 'vitro'
pcs$kind <- 'MA'
pcs[grepl('.M', as.character(pcs$sample), fixed=T),]$env <- 'vivo'
pcs[grepl('LM', as.character(pcs$sample), fixed=T),]$kind <- 'bulk'
pcs$tt <- 'T1'
pcs[pcs$time=='01',]$tt <- 'T0'
pcs[pcs$time=='0',]$tt <- 'basale'
table(pcs$tt)

ggplot(data = pcs, aes(x = PC1, y = PC2, color = time)) +geom_point(size = 4) + xlab(paste0("PC1", ': ', round(percentVar[1] *100), "% variance")) + ylab(paste0("PC2", ":", round(percentVar[2] *100), "% variance")) +theme_bw()

ff <- function(pc1, pc2) {
pcs <- data.frame(PC1 = pca$x[, pc1], PC2 = pca$x[, pc2], sample = rownames(pca$x))
pcs$model <- (unlist(lapply(strsplit(as.character(pcs$sample),'.', fixed=TRUE), function(x){ x[1] })))
pcs$model2 <- substr(pcs$model, 0,7)
pcs$time <- unlist(lapply(strsplit(as.character(pcs$sample),'.', fixed=T), function(x){ if (length(x)>2) {x[3]} else {'B'} }))
ggplot(data = pcs, aes_string(x = 'PC1', y = 'PC2', color='model2', shape = 'time')) +geom_point(size = 4) + xlab(paste0("PC", pc1,': ', round(percentVar[pc1] *100), "% variance")) + ylab(paste0("PC", pc2, ":", round(percentVar[pc2] *100), "% variance")) +theme_bw()+scale_color_manual(values=cbPalette2[-1])

}