# Heterogeneity and evolution of DNA mutation rates in microsatellite-stable colorectal cancer

Snakemake rules to reproduce all main Figures panels (with the exception of 1A, made on biorender) can be found in the directory `dataset_Figures_Tables`.
Keep in mind that final figures were assembled manually (e.g. for legend positioning) and that some rules produce more than one plot (with/without legend) to ease
those manual fixes.

List of panels:
- fig_1b_MR.svg                                                                                                                                                                  
- fig_1c_dnds.svg
- fic_2cbis_cosmic.pdf                                                                                                                                                               
- fig_2a_cosmic.pdf                                                                                                                                                                
- fig_2abis_cosmic.pdf                                                                                                                                                             
- fig_2b_cosmic.pdf                                                                                                                                                                    
- fig_2c_cosmic.pdf                                                                                                                                                                
- fic_2cbis_cosmic.pdf
- fig_3a_MR.svg                                                                                                                                                                   
- fig_3b_wilcox_mb.svg                                                                                                                                                            
- fig_4b_CN.svg                                                                                                                                                                    
- fig_4c_cor.svg                                                                                                                                                              
- fig_5a_subclonal.svg                                                                                                                                                                                                                                    
- fig_5b_cum_nolegend.svg                                                                                                                                                      
- fig_5c_slopes_pairedscatter.svg                                                                                                                                                   

Also Extended Data Figures and Tables can be reproduced in the directory `dataset_Figures_Tables`. Keep in mind that alternative versions were generated for some figures and rules producing them are still available.
The numbering of Extended materials does not always reflect the one in the manuscript, but these exceptions are commented in the Snakefile and fill be fixed for the final revised manuscript (^_^).

A Dockerfile with the same R version, R packages and overall linux environment used to produce the Figures is available here: https://github.com/vodkatad/snakemake_docker/blob/master/Dockerfiles/godot/Dockerfile 
