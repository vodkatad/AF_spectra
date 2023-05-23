
setwd('/home/data/Dropbox/work/evol/MA/sign_input_EDT5/')
library(signature.tools.lib)
vcf_files <- c("CRC0282LMO-0-B.pass.vcf.gz","CRC0282_vitroMA.vcf.gz","CRC0282_vivoMA.vcf.gz","CRC0327LMO-0-B.pass.vcf.gz","CRC0327_vitroMA.vcf.gz","CRC0327_vivoMA.vcf.gz","CRC0441LMO-0-B.pass.vcf.gz","CRC0441_vitroMA.vcf.gz","CRC0441_vivoMA.vcf.gz","CRC1078LMO-0-B.pass.vcf.gz","CRC1078_vitroMA.vcf.gz","CRC1078_vivoMA.vcf.gz","CRC1307LMO-0-B.pass.vcf.gz","CRC1307_vitroMA.vcf.gz","CRC1307_vivoMA.vcf.gz","CRC1502LMO-0-B.pass.vcf.gz","CRC1502_vitroMA.vcf.gz","CRC1502_vivoMA.vcf.gz","CRC1599LMO-0-B.pass.vcf.gz","CRC1599LM_vitroMA.vcf.gz","CRC1599LM_vivoMA.vcf.gz","CRC1599PRO-0-B.pass.vcf.gz","CRC1599PR_vitroMA.vcf.gz")

sample_names <- gsub('.vcf.gz', '', vcf_files)
names(vcf_files) <- sample_names

get_catalogue <- function(file, name, genome='hg38') {
  res <- vcfToSNVcatalogue(file, genome.v=genome)
  cat <- res$catalogue
  colnames(cat) <- name
  return(cat)
}

cat <- mapply(get_catalogue, vcf_files, sample_names)
SNV_catalogues <- do.call(cbind,cat)

plotSubsSignatures(signature_data_matrix = SNV_catalogues, plot_sum = TRUE, output_file='subs_sign.pdf')

fit <- Fit(catalogues = SNV_catalogues,
    signatures = getOrganSignatures("Colorectal"),
    useBootstrap = TRUE,
    nboot = 100,
    nparallel = 4)
plotFit(fit,outdir = "signatureFitTest/")
snv_exp <- fit$exposures

comment out all bastards:
if(tools:::.BioC_version_associated_with_R_version()<3.5){
    gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms))
  }else{
    gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms),pruning.mode = "coarse")
  }
  
  we are < 3.5 ?
  
but if we keep this if branch:
snv_exp <- fit$exposures
Error in GenomeInfoDb:::getDanglingSeqlevels(x, new2old = new2old, pruning.mode = pruning.mode,  : 
  The following seqlevels are to be dropped but are currently in use
  (i.e. have ranges on them): chrM, chr1_GL383518v1_alt,
  chr1_GL383519v1_alt, chr1_GL383520v2_alt, chr1_KI270759v1_alt,
  chr1_KI270760v1_alt, chr1_KI270761v1_alt, chr1_KI270762v1_alt,
  chr1_KI270763v1_alt, chr1_KI270764v1_alt, chr1_KI270765v1_alt,
  chr1_KI270766v1_alt, chr1_KI270892v1_alt, chr2_GL383521v1_alt,
  chr2_GL383522v1_alt, chr2_GL582966v2_alt, chr2_KI270767v1_alt,
  chr2_KI270768v1_alt, chr2_KI270769v1_alt, chr2_KI270770v1_alt,
  chr2_KI270771v1_alt, chr2_KI270772v1_alt, chr2_KI270773v1_alt,
  chr2_KI270774v1_alt, chr2_KI270775v1_alt, chr2_KI270776v1_alt,
  chr2_KI270893v1_alt, chr2_KI270894v1_alt, chr3_GL383526v1_alt,
  chr3_JH636055v2_alt, chr3_KI270777v1_alt, chr3_KI270778v1_alt,
  chr3_KI270779v1_alt, chr3_KI270780v1_alt, chr3_KI270781v1_alt,
  chr3_KI270782v1_alt, chr3_KI270783v1_alt, chr3_KI270784v1_alt,
  chr3_KI270895v1_alt, chr3_KI270924v1_alt, chr3_KI270934v1_alt,
  chr3_KI270935v1_alt,


kept the other.

fit <- Fit(catalogues = SNV_catalogues,
    signatures =  getSignaturesForFitting(
       "Colorectal",
       typemut = "subs",
       commontier = "T2"),
    useBootstrap = TRUE,
    nboot = 100,
    nparallel = 4)
plotFit(fit,outdir = "signatureFitTestT2/")

Error in { : 
  task 1 failed - "'list' object cannot be coerced to type 'double'"


pORCODIOMAIALE

others: sigflow (very standard), sigfit (very different since bayesian), deconstructSigs (famous but hg19?)


...
options(mc.cores = 10)
library(sigfit)
data("cosmic_signatures_v2")
count <- t((SNV_catalogues))

mcmc_samples_fit <- fit_signatures(counts = count, 
                                   signatures = cosmic_signatures_v2,
                                   iter = 2000, 
                                   warmup = 1000, 
                                   chains = 1, 
                                   seed = 1756)
                                   
                                   Fitting 30 signatures using multinomial model
---
Stan sampling:
SAMPLING FOR MODEL 'sigfit_fit' NOW (CHAIN 1).
Chain 1: 
Chain 1: Gradient evaluation took 0.000792 seconds
Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 7.92 seconds.
Chain 1: Adjust your expectations accordingly!
Chain 1: 
Chain 1: 
Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
Chain 1: 
Chain 1:  Elapsed Time: 40.422 seconds (Warm-up)
Chain 1:                15.009 seconds (Sampling)
Chain 1:                55.431 seconds (Total)
Chain 1: 
Warning message:
Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
Running the chains for more iterations may help. See
https://mc-stan.org/misc/warnings.html#tail-ess 

exposures <- retrieve_pars(mcmc_samples_fit, 
                           par = "exposures", 
                           hpd_prob = 0.90)
                           
write.table(exposures, file="sigfit.tsv", sep="\t", quote=F)

