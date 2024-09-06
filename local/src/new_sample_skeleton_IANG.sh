#!/bin/bash

mkdir $1
#mkdir $1/platypus_nobin_00/
mkdir $1/platypus_nobin_indels_00 $1/sequenza $1/mutect_nobin $1/platypus_nobin_00
#mkdir $1/platypus_nobin_indels_00 
#mkdir -p $1/mutect_VAF
#mkdir -p $1/univMutect
#cd $1/platypus_nobin_00
cd $1/platypus_nobin_00/
if [ "$1" == "CRC1599PR" ]
then
    ln -s ../../../local/share/snakerule/Snakefile_clones_platypus_polished_nobinomial_00_bugbedtools Snakefile
else
    ln -s ../../../local/share/snakerule/Snakefile_clones_platypus_polished_nobinomial_00_notopi Snakefile
fi
ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk

#cd $1/platypus_nobin_indels_00/
cd ../mutect_nobin/
ln -s ../../../local/share/snakerule/Snakefile_clones_real_segments Snakefile
ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk
#cd $1/mutect_VAF/
#ln -s ../../../local/share/snakerule/Snakefile_clones_real_segments_VAFthr Snakefile
#ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk
#cd $1/univMutect
#ln -s ../../../local/share/snakerule/Snakefile_univMutect_pileup Snakefile
#ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk

cd ../sequenza
ln -s ../../../local/share/snakerule/Snakefile_sequenza Snakefile
ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk

#cd ../MutationalPatterns/
#ln -s ../../../local/share/snakerule/Snakefile_mutationalpatterns Snakefile
#ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk

#mkdir $1/tree
#cd $1/tree
#ln -s ../../../local/share/snakerule/Snakefile_tree_v2 Snakefile
#ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk
#mkdir $1/sequenza
#cd $1/sequenza
#ln -s ../../../local/share/snakerule/Snakefile_sequenza Snakefile
#ln -s ../../../local/share/snakemake/conf_$1_mutect_platypus_polished.sk conf.sk