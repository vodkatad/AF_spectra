# Heterogeneity and evolution of DNA mutation rates in microsatellite-stable colorectal cancer

Snakemake rules to reproduce all main Figures panels can be found in the directory `dataset_Figures_Tables`.
Keep in mind that final assembled figures were manually fixed (e.g. for legend positioning) and that some rules produces more than one plot (with/without legend) to ease
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
- fig_5a_subclonal.svg                                                                                                                                                                                                                                                                                                                      100%   16KB 734.9KB/s   00:00    
- fig_5b_cum_nolegend.svg                                                                                                                                                      
- fig_5c_slopes_pairedscatter.svg                                                                                                                                                   

Also Extended Data Figures and Tables can be obtained in the directory `dataset_Figures_Tables`. Keep in mind that some alternative versions for the figures exists and rules producing them are still available.
The numbering of Extended Material does not always reflect the one in the manuscript, but these exceptions were commented in the Snakefile.
A Dockerfile with the same R version and R packages used to produce the Figures is available here: https://github.com/vodkatad/snakemake_docker/blob/master/Dockerfiles/godot/Dockerfile 
