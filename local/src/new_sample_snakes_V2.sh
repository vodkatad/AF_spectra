#!/bin/bash
#shopt -s expand_aliases
#alias lsnakemake='/home/egrassi/.local/bin/snakemake --log-handler-script /home/egrassi/sysadm/snakemake_slack.py'

#cd $1/platypus_nobin_00
#snakemake -j 12 all.MR_ov all_gain_vcf all_gained_named.tsv 
#snakemake -j 12 dnds.tsv dndsvitro.tsv
#snakemake -j 12 all_mrca.tsv
#snakemake vitro.merged.vcf.gz vivo.merged.vcf.gz

#cd $1/platypus_nobin_indels_00/
#cd ../platypus_nobin_indels_00/
#rm *tsv *gz *ov *png *log *multiallelic
#snakemake -j 12 all.MR_ov

#cd ../mutect_nobin/
#snakemake all_R

#cd ../sequenza
#snakemake -j 12 merged_heatmap.png

#cd ../MutationalPatterns/
#snakemake Homo_sapiens.GRCh37.75_autosomal_exon_merged_sorted.enrich.tsv


cd $1/tree/
snakemake tree_bulk_vitro.pdf nobs_tree_bulk_vitro.svg