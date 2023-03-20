# Libraries
if (!require(vcfR)) stop("Package 'vcfR' missing\n.")

vt = 0.1
bs = 10000
#bs = 1

# Read in the vcf
vcf = read.vcfR('/scratch/trcanmed/AF_spectra/local/share/data/third_shipment_bulk/CRC1599LMO-0-B.pass.vcf.gz')
afLM     = extract.gt(vcf, element = "AF", as.numeric = TRUE)

vcf = read.vcfR('/scratch/trcanmed/AF_spectra/local/share/data/third_shipment_bulk/CRC1599PRO-0-B.pass.vcf.gz')
# What is the depth?
afPR     = extract.gt(vcf, element = "AF", as.numeric = TRUE)

exsubclM <- afLM[afLM[,1] > 0.05 & afLM[,1] < 0.24, 1]
hist(exsubclM, breaks=50, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col='#ff9900', main="", xlim=c(0, 0.24), ylim=c(0, 500))


exsubclP <- afPR[afPR[,2] > 0.05 & afPR[,2] < 0.24, 1]
hist(exsubclP, breaks=25, cex=1.5, xlab="Allelic frequency (f)", ylab="Number of muts", border="black", col='#ffff00', main="", xlim=c(0, 0.24))
