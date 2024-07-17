library(ggplot2)
d <- data.frame(x=runif(100))

ggplot(data=d, aes(x=x))+geom_density()
ggplot(data=d, aes(x=x, y=after_stat(count)/nrow(d)))+geom_density()+geom_histogram()
ggplot(data=d, aes(x=x, y=after_stat(count)))+geom_density()+geom_histogram()
ggplot(data=d, aes(x=x, y=after_stat(scaled)))+geom_density()

ggplot(data=d, aes(x=x, y=after_stat(density)))+geom_density()
