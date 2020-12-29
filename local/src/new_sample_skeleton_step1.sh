#!/bin/bash

mkdir $1
mkdir $1/platypus_nobin_indels $1/sequenza $1/mutect_nobin $1/platypus_nobin $1/MutationalPatterns
cd $1/platypus_nobin
ln -s ../../../local/share/snakerule/Snakefile_clones_platypus_polished_nobinomial Snakefile
cp ../../../local/share/snakemake/conf_CRC1502_clones_all_mutect_platypus_polished.sk ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk
ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk
cd ../platypus_nobin_indels/
ln -s ../../../local/share/snakerule/Snakefile_clones_platypus_polished_nobinomial Snakefile
ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished_indels.sk conf.sk
cd ../mutect_nobin/
ln -s ../../../local/share/snakerule/Snakefile_clones_real_segments Snakefile
ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk
cd ../sequenza
ln -s ../../../local/share/snakerule/Snakefile_sequenza Snakefile
ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk
cd ../MutationalPatterns/
ln -s ../../../local/share/snakerule/Snakefile_mutationalpatterns Snakefile
ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk

echo "Now add your samples info in ../local/share/snakemake/conf_$1_mutect_platypus_polished.sk, create ../local/share/snakemake/conf_$1_mutect_platypus_polished_indels.sk for indels - then call ../local/src/new_sample_snakes.sh"
