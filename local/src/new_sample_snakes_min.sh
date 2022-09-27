#!/bin/bash
#shopt -s expand_aliases

cd $1/platypus_nobin
snakemake -j 12 all.MR_ov

cd ../platypus_nobin_indels/
snakemake -j 12 all.MR_ov

