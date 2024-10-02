library(ggplot2)

set.seed(42)
flat <- runif(300)
low <- runif(50, min=0, max=0.2)
high <- runif(300, min=0.8, max=1)
flatless <- runif(100)

d <- data.frame(x=c(flat, low, high, flatless), type=c(rep('flat', 300), rep('low', 50), rep('high', 300), rep('flatless', 100)))

ggplot(data=d, aes(x=x, fill=type))+geom_histogram(position='identity', alpha=0.4)
ggplot(data=d, aes(x=x, color=type))+geom_density()
ggplot(data=d, aes(x=x, color=type))+geom_density(aes(y=after_stat(count)))
ggplot(data=d, aes(x=x, color=type))+geom_density(aes(y=after_stat(scaled)))

ggplot(data=d, aes(x=x, color=type))+geom_density(aes(y=after_stat(scaled*n)))
ggplot(data=d, aes(x=x, color=type))+geom_density(aes(y=after_stat(count*1/30)))

ggplot(data=d, aes(x=x, color=type))+geom_density(aes(y=after_stat(count*0.1)))
ggplot(data=d, aes(x=x, fill=type))+geom_histogram(position='identity', alpha=0.4, binwidth=0.1)
ggplot(data=d, aes(x=x, color=type))+geom_density(aes(y=after_stat(count*0.1)))+geom_histogram(position='identity', alpha=0.4, binwidth=0.1)


ggplot(data=d, aes(x=x, color=type))+geom_density(aes(y=after_stat(count*0.01)))
ggplot(data=d, aes(x=x, fill=type))+geom_histogram(position='identity', alpha=0.4, binwidth=0.01)

#     after_stat(density)
#density estimate.

#after_stat(count)
#density * number of points - useful for stacked density plots.

#after_stat(scaled)
#density estimate, scaled to maximum of 1.

#after_stat(n)
#number of points.

#after_stat(ndensity)
#alias for scaled, to mirror the syntax of stat_bin().

#https://stackoverflow.com/questions/65631194/scale-density-curve-made-with-geom-density-to-similar-height-of-geom-histogram
