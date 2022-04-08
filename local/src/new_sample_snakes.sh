#!/bin/bash
#shopt -s expand_aliases
#alias lsnakemake='/home/egrassi/.local/bin/snakemake --log-handler-script /home/egrassi/sysadm/snakemake_slack.py'

cd $1/platypus_nobin
snakemake -j 12 all.MR_ov all_gain_vcf all_gained_named.tsv 
snakemake -j 12 dnds.tsv dndsvitro.tsv

cd ../platypus_nobin_indels/
snakemake -j 12 all.MR_ov

cd ../mutect_nobin/
snakemake all_R

cd ../sequenza
snakemake -j 12 merged_heatmap.png

cd ../MutationalPatterns/
snakemake Homo_sapiens.GRCh37.75_autosomal_exon_merged_sorted.enrich.tsv
