setwd('/mnt/trcanmed/snaketree/prj/AF_spectra/dataset')
d <- read.table('all_bexp_different_coverage.tsv', header=T, sep="\t")
dd <- d[d$coverage=="10x",]

sum4 <- sapply(unique(dd$sample), function(x) { sum(dd[dd$generation==4 & dd$sample==x,'called_in_gen'])})
sum3 <- sapply(unique(dd$sample), function(x) { sum(dd[dd$generation==3 & dd$sample==x,'called_in_gen'])})
sum2 <- sapply(unique(dd$sample), function(x) { sum(dd[dd$generation==2 & dd$sample==x,'called_in_gen'])})
sum1 <- sapply(unique(dd$sample), function(x) { sum(dd[dd$generation==1 & dd$sample==x,'called_in_gen'])})
calls <- data.frame(gen1=sum1, gen2=sum2, gen3=sum3, gen4=sum4)

mm <- melt(calls)
mm$clone <- unique(dd$sample)
ggplot(mm, aes(x=variable, y=value, group=clone, color=clone)) +
geom_point()+geom_line()

dd <- d[d$coverage=="20x",]

sum4 <- sapply(unique(dd$sample), function(x) { sum(dd[dd$generation==4 & dd$sample==x,'called_in_gen'])})
sum3 <- sapply(unique(dd$sample), function(x) { sum(dd[dd$generation==3 & dd$sample==x,'called_in_gen'])})
sum2 <- sapply(unique(dd$sample), function(x) { sum(dd[dd$generation==2 & dd$sample==x,'called_in_gen'])})
sum1 <- sapply(unique(dd$sample), function(x) { sum(dd[dd$generation==1 & dd$sample==x,'called_in_gen'])})
calls <- data.frame(gen1=sum1, gen2=sum2, gen3=sum3, gen4=sum4)

mm <- melt(calls)
mm$clone <- unique(dd$sample)
ggplot(mm, aes(x=variable, y=value, group=clone, color=clone)) +
  geom_point()+geom_line()