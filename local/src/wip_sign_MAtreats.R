load('/scratch/trcanmed/AF_spectra/dataset_MAtreats/sign.Rdata')
library(MutationalPatterns)
library('BSgenome.Hsapiens.UCSC.hg38')
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)


nt_type_occurrences <- mut_type_occurrences(vcfs[3], ref_genome)
t_type_occurrences <- mut_type_occurrences(vcfs[5], ref_genome)
plot_spectrum(t_type_occurrences)
plot_spectrum(nt_type_occurrences)
plot_96_profile(mut_mat)
chromosomes <- seqnames(get(ref_genome))[1:22]

plot_rainfall(vcfs[[5]],
              title = names(vcfs[5]),
              chromosomes = chromosomes, cex = 1.5, ylim = 1e+09
)

plot_rainfall(vcfs[[3]],
              title = names(vcfs[3]),
              chromosomes = chromosomes, cex = 1.5, ylim = 1e+09
)

# strand bias: nope!


genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
strand <- mut_strand(vcfs[[5]], genes_hg38)
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes_hg38)

plot_192_profile(mut_mat_s[, c(3,5)])
strand_counts <- strand_occurrences(mut_mat_s)
strand_bias <- strand_bias_test(strand_counts)

# replicative stress


# regions

