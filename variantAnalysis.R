library(vcfR)
library(adegenet)
library("SNPRelate")

# Read vcf file
vcf <- read.vcfR("C:\\Users\\Javan\\Desktop\\sars-cov-2Project\\variantCalling\\SRX8540989.filtered.vcf")

#Shows give the vcf report
vcf

#Converting VCF data to a genlight object
x <- vcfR2genlight(vcf)
x

# Check the genotype
# vcfR
gt <- extract.gt(vcf, element = "GT")

t = t(as.matrix(x))

View(t)

#PLOTTING PCA
vcf.fn<-"C:\\Users\\Javan\\Desktop\\sars-cov-2Project\\variantCalling\\SRX8540989.filtered.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2] ,col=as.numeric(substr(ccm_pca$sample, 1,3) == 'CCM')+3, pch=2)

jColors <- c('chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3',
             'mediumorchid2', 'turquoise3', 'wheat4', 'slategray2') ## choosing specific colors for your plot!

plot(tab_popr$EV2, tab_popr$EV1, col=jColors[as.integer(tab_popr$pop)], xlab = "EV2", ylab = "EV1")
legend("topleft", legend = levels(tab_popr$pop), pch = "o", col = jColors[1:nlevels(tab_popr$pop)])

