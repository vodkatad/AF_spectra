#!/bin/bash
#shopt -s expand_aliases
#alias lsnakemake='/home/egrassi/.local/bin/snakemake --log-handler-script /home/egrassi/sysadm/snakemake_slack.py'

cd $1/platypus_nobin_00

#rm *_ov *ov.tsv
#snakemake -j 12 all.MR_ov
#snakemake -j 12 all.MR_ov all_gain_vcf all_gained_named.tsv 
#snakemake -j 12 dnds.tsv dndsvitro.tsv
#snakemake -j 12 all_mrca.tsv
#snakemake vitro.merged.vcf.gz vivo.merged.vcf.gz
#snakemake repliseq_gained.genomicregions.svg
#snakemake vivo.merged.vcf.gz
#snakemake mutinfo.tsv.gz
snakemake  all.length
#cd ../platypus_nobin_indels_00/
#rm *_ov *ov.tsv
#snakemake -j 12 all.MR_ov
#cd ../platypus_nobin_indels_00/
#rm *tsv *gz *ov *png *log *multiallelic
#snakemake -j 12 all.MR_ov
#snakemake dnds_double_vitro.tsv
#snakemake mutinfo.tsv.gz
#cd ../mutect_nobin/
#snakemake all_mrca.tsv
#snakemake all.MR_ov
#snakemake all_R
#cd $1/mutect_VAF/
#snakemake -j 6 all.MR_ov
#cd $1/univMutect/
#snakemake -j 6 all.MR_ov
#snakemake -j 12 all_mrca.tsv
#cd $1/mutect_nobin
#snakemake subclonal_mr_corr.png
#cd ../sequenza
#snakemake -j 12 merged_heatmap.png
#cd $1/mutect_nobin/
#snakemake -j 12 all_T0_fbcalls.tsv.gz

#cd ../MutationalPatterns/
#snakemake Homo_sapiens.GRCh37.75_autosomal_exon_merged_sorted.enrich.tsv


#cd $1/tree/
#snakemake tree_bulk_vitro.pdf nobs_tree_bulk_vitro.svg
#snakemake bulkmutinfo.tsv.gz

#cd $1/tree
#snakemake -f bulk.nonsyn.binary.tsv.gz