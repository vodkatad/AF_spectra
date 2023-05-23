bulk_f  <- snakemake@input[['bulk']]
gained_f  <- snakemake@input[['gained']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape)
load(theme)

d <- read.table(bulk_f, header=TRUE, fill=TRUE)
d_bulk <- as.data.frame(t(data.frame(row.names=colnames(d), n=d[,1])))

d <- read.table(gained_f, header=TRUE, fill=TRUE)
d_gained <- as.data.frame(t(data.frame(row.names=colnames(d), n=d[,1])))

ess <- c(d_gained$n_gainedess, d_bulk$n_gainedess)
total <- c(d_gained$n_gained, d_bulk$n_gained)

ns <- data.frame(row.names=c('de_novo', 'truncal'), essential=ess,
                 other=total-ess)

sink(log_f)
print('chisq')
fisher.test(ns)
sink()

save.image(paste0(outplot, '.Rdata'))

ns$id <- rownames(ns)
pd <- melt(ns)
df <- pd %>% 
  group_by(id) %>% # Variable to be transformed
  mutate(perc = `value` / sum(`value`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

sink(log_f, append=TRUE)
print(df)
sink()


#ggplot(df, aes(x = "", y = perc, fill = annot.type)) +
#  geom_col() +
#  coord_polar(theta = "y")+facet_grid(~SNV)+theme_light()

pal <- brewer.pal(n = 8, name = "Dark2")

colnames(df) <- c('SNV', 'gene', 'n', 'perc', 'label')
p <- ggplot(df, aes(x = "", y = perc, fill = gene)) +
  geom_col() + facet_grid(~SNV)+ unmute_theme+scale_fill_manual(values=pal[c(4,5)]) + scale_y_continuous(expand = c(0, 0))

ggsave(outplot, plot=p, width=89, height=56, units="mm")
