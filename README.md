# Heterogeneity and evolution of DNA mutation rates in microsatellite-stable colorectal cancer

Snakemake rules to reproduce all main Figures panels (with the exception of 1A, made on biorender) can be found in the directory `dataset_Figures_Tables`.
Keep in mind that final figures were assembled manually (e.g. for legend positioning) and that some rules produce more than one plot (with/without legend) to ease
those manual fixes.                                                                                                                                           

Also Extended Data Figures and Tables can be reproduced in the directory `dataset_Figures_Tables`. Keep in mind that alternative versions were generated for some figures and rules producing them are still available.
The numbering of Extended materials does not always reflect the one in the manuscript, but these exceptions are commented in the Snakefile and fill be fixed for the final revised manuscript (^_^).

A Dockerfile with the same R version, R packages and overall linux environment used to produce the Figures is available here: https://github.com/vodkatad/snakemake_docker/blob/master/Dockerfiles/godot/Dockerfile and on the docker hub (https://hub.docker.com/r/egrassi/godot):
```
$ docker pull docker.io
$ git clone git@github.com:vodkatad/AF_spectra.git
$ cd AF_spectra
$ WD=$PWD
$ cd dataset_Figures_Tables/
$ docker run -v $WD:$WD -w $PWD -it egrassi/godot  
# snakemake fig_1b_MR.svg
```
