library(ggplot2)
library(reshape)
data <- read.table('/scratch/trcanmed/AF_spectra/local/share/data/alltogether_clocks.tsv', sep="\t", header=T, stringsAsFactors = FALSE)

colnames(data) <- c('name', '1e-4', '5e-5', '1e-5', '5e-6', 'clone')
data <- data[, c(1, 6, 5,4,3,2)]
us <- data.frame(clone=c('clone9', 'clone8', 'clone2'), b=c(1.2229,1.2686,1.1886), d=c(0.7312, 0.7385, 0.7346),
                 gr=c(0.4916,0.5300,0.4540))


data$clock <- 'old'
data[data$clone=="55C-9_CHAT", 'clock']  <- "CHAT"
data[data$clone=="55C-9_CHAT", 'clone']  <- "clone9"
data[data$clone=="55C-9_TFAP2B", 'clock']  <- "TFAP2B"
data[data$clone=="55C-9_TFAP2B", 'clone']  <- "clone9"
data[data$clone=="58C-9_IRX2", 'clock']  <- "IRX2"
data[data$clone=="58C-9_IRX2", 'clone']  <- "clone9"
data[data$clone=="58C-9_MYSTERY", 'clock']  <- "MYSTERY"
data[data$clone=="58C-9_MYSTERY", 'clone']  <- "clone9"

data[data$clone=="55C-9c_CHAT", 'clock']  <- "CHAT"
data[data$clone=="55C-9c_CHAT", 'clone']  <- "clone9c"
data[data$clone=="55C-9c_TFAP2B", 'clock']  <- "TFAP2B"
data[data$clone=="55C-9c_TFAP2B", 'clone']  <- "clone9c"
data[data$clone=="58C-9c_IRX2", 'clock']  <- "IRX2"
data[data$clone=="58C-9c_IRX2", 'clone']  <- "clone9c"
data[data$clone=="58C-9c_MYSTERY", 'clock']  <- "MYSTERY"
data[data$clone=="58C-9c_MYSTERY", 'clone']  <- "clone9c"

birth <- data[data$name=="birth rate",]
birth$name <- NULL
birthlong <- melt(birth, id.vars=c("clone", "clock"))
#birthlong$variable <- as.numeric(as.character(birthlong$variable))
birthlong <- birthlong[birthlong$clone %in% c('clone9', 'clone9c'),]
ggplot(data=birthlong, aes(y=value, color=clone, x=variable))+geom_point()+
  theme_bw()+facet_grid(~clock)+geom_hline(yintercept=us[us$clone=="clone9", 'b'])
  

birth <- data[data$name=="death rate",]
birth$name <- NULL
birthlong <- melt(birth, id.vars=c("clone", "clock"))
birthlong <- birthlong[birthlong$clone %in% c('clone9', 'clone9c'),]
ggplot(data=birthlong, aes(y=value, color=clone, x=variable))+geom_point()+
  theme_bw()+facet_grid(~clock)+geom_hline(yintercept=us[us$clone=="clone9", 'd'])

birth <- data[data$name=="birth rate",]
birth$name <- NULL
birthlong <- melt(birth, id.vars=c("clone", "clock"))
m <- merge(birthlong, us, by="clone")
ggplot(data=m, aes(y=value, color=clone, x=variable))+geom_point()+
     theme_bw()+facet_grid(~clock)+geom_hline(aes(yintercept=b, color=clone))


birth <- data[data$name=="death rate",]
birth$name <- NULL
birthlong <- melt(birth, id.vars=c("clone", "clock"))
m <- merge(birthlong, us, by="clone")
ggplot(data=m, aes(y=value, color=clone, x=variable))+geom_point()+
  theme_bw()+facet_grid(~clock)+geom_hline(aes(yintercept=d, color=clone))


birth <- data[data$name=="net growhthrate",]
birth$name <- NULL
birthlong <- melt(birth, id.vars=c("clone", "clock"))
m <- merge(birthlong, us, by="clone")
ggplot(data=m, aes(y=value, color=clone, x=variable))+geom_point()+
  theme_bw()+facet_grid(~clock)+geom_hline(aes(yintercept=gr, color=clone))


