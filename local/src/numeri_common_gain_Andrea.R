library(ggplot2)
library(reshape)
d <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/stat_gain_common_loss_noall', sep="\t", stringsAsFactors = FALSE)
colnames(d) <- c('id', 'n','class')

d$clone <- sapply(d$id, function(x) {y<-strsplit(x, '_')[[1]][3]; return(y[1])})
d$clone <- gsub('00/','', d$clone, fixed=TRUE)

n_common_gain <- function(filter_data, data, filter_column, remove) {
  mdata <- data[data$class != remove,]
  mdata <- mdata[mdata[,filter_column] == filter_data,]
  tot <- sum(mdata$n)
  res <- c()
  for (cl in unique(mdata$class)) {
    res <- c(res, mdata[mdata$class== cl, 'n'] / tot)
    
  }
  names(res) <- unique(mdata$class)
  return(res)
}

# cagate ma per le percentuali serve
ll <- t(sapply(unique(d$clone), n_common_gain, d, 'clone', 'loss'))
ll <- as.data.frame(ll)

ll$id <- rownames(ll)

mpd <- melt(ll, id='id')
mpd$variable <- as.factor(mpd$variable, levels=c('common', 'gained'))
ggplot(data=mpd, aes(x=id, y=value, fill=variable))+geom_col(position='stack')+theme_bw(base_size=20)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+scale_fill_manual(values=c('darkgreen', 'darkgoldenrod'), name="Mutation type")+ylab('N.')

## con loss

ll <- t(sapply(unique(d$clone), n_common_gain, d, 'clone', 'NOPE'))
ll <- as.data.frame(ll)

ll$id <- rownames(ll)

mpd <- melt(ll, id='id')
mpd$variable <- as.factor(mpd$variable, levels=c('common', 'gained', 'loss'))
ggplot(data=mpd, aes(x=id, y=value, fill=variable))+geom_col(position='stack')+theme_bw(base_size=20)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+scale_fill_manual(values=c('darkgreen', 'darkgoldenrod', 'darkred'), name="Mutation type")+ylab('N.')


#### per modello
d$model <- sapply(d$clone, function(x) {y<-strsplit(x, '-')[[1]][1]; return(y[1])})

sum_n_common_gain <- function(filter_data, data, filter_column, remove) {
  mdata <- data[data$class != remove,]
  mdata <- mdata[mdata[,filter_column] == filter_data,]
  tot <- sum(mdata$n)
  res <- c()
  for (cl in unique(mdata$class)) {
    mysum <- sum(mdata[mdata$class== cl, 'n'])
    res <- c(res, mysum / tot)
    
  }
  names(res) <- unique(mdata$class)
  return(res)
}


# cagate ma per le percentuali serve
ll <- t(sapply(unique(d$model), sum_n_common_gain, d, 'model', 'loss'))
ll <- as.data.frame(ll)

ll$id <- rownames(ll)

mpd <- melt(ll, id='id')
mpd$variable <- as.factor(mpd$variable, levels=c('common', 'gained'))
ggplot(data=mpd, aes(x=id, y=value, fill=variable))+geom_col(position='stack')+theme_bw(base_size=20)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+scale_fill_manual(values=c('darkgreen', 'darkgoldenrod'), name="Mutation type")+ylab('N.')


ll <- t(sapply(unique(d$model), sum_n_common_gain, d, 'model', 'NOPE'))
ll <- as.data.frame(ll)

ll$id <- rownames(ll)

mpd <- melt(ll, id='id')
mpd$variable <- as.factor(mpd$variable, levels=c('common', 'gained', 'loss'))
ggplot(data=mpd, aes(x=id, y=value, fill=variable))+geom_col(position='stack')+theme_bw(base_size=20)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+scale_fill_manual(values=c('darkgreen', 'darkgoldenrod', 'darkred'), name="Mutation type")+ylab('N.')





