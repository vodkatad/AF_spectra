data_f  <- snakemake@input[['data']]

log_f <- snakemake@log[['log']]
outplot <- snakemake@output[['plot']]
theme <- snakemake@input[['theme']]
save.image(paste0(outplot, '.Rdata'))

library(dplyr)
library(ggplot2)
load(theme)

pd <- read.table(data_f, sep="\t", stringsAsFactors=FALSE, header=TRUE)

ns <- data.frame(row.names=pd[pd$SNV == "de_novo", "annot.type"], gained=pd[pd$SNV == "de_novo", "n"],
                 bulk=pd[pd$SNV == "truncal", "n"])


sink(log_f)
print('chisq')
chisq.test(ns)
sink()

save.image(paste0(outplot, '.Rdata'))

df <- pd %>% 
  group_by(SNV) %>% # Variable to be transformed
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

sink(log_f, append=TRUE)
print(df)
sink()


#ggplot(df, aes(x = "", y = perc, fill = annot.type)) +
#  geom_col() +
#  coord_polar(theta = "y")+facet_grid(~SNV)+theme_light()

p <- ggplot(df, aes(x = "", y = perc, fill = annot.type)) +
  geom_col() + facet_grid(~SNV)+ unmute_theme+scale_fill_brewer(palette="Dark2")

ggsave(outplot, plot=p, width=89, height=56, units="mm")
