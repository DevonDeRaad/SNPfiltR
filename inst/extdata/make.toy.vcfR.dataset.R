devtools::install_github("DevonDeRaad/SNPfiltR")

#clean package
#library(pkgdown)
#clean_site(pkg = "SNPfiltR")
#library(devtools)
#clean_vignettes(pkg = "SNPfiltR")

#load packages
library(vcfR)
library(SNPfiltR)

#read in vcf and make subset example datasets
vcffog<-read.vcfR("~/Desktop/hipposideros/n3.vcf")

#subset to 500 random SNPs
vcfR.example<-vcffog[sample.int(151015, 500), c(1:20)]

#save this vcf
write.vcf(vcfR.example, file = "~/Downloads/SNPfilter.example.vcf.gz")

#save this vcfR as available data for the SNPfiltR package
usethis::use_data(vcfR.example, internal = TRUE, overwrite = T)

#clean global environment
rm(list=ls())

#load and test examples

#missing by snp
missing.by.snp(vcfR = SNPfiltR:::vcfR.example)
missing.by.snp(vcfR = SNPfiltR:::vcfR.example, cutoff = .6)

#missing by sample
missing.by.sample(vcfR = SNPfiltR:::vcfR.example)
missing.by.sample(vcfR = SNPfiltR:::vcfR.example, cutoff = .7)

#hard filter
hard.filter.vcf(vcfR = SNPfiltR:::vcfR.example, depth = 5)
hard.filter.vcf(vcfR = SNPfiltR:::vcfR.example, depth = 5, gq = 30)

#max depth
max_depth(vcfR = SNPfiltR:::vcfR.example)
max_depth(vcfR = SNPfiltR:::vcfR.example, maxdepth = 100)

#allele balance
filter.allele.balance(vcfR = SNPfiltR:::vcfR.example)

#distance thin
#this function needs work
#must build in a fail if no distance is specified
#must fix the issue that de novo assembled loci are not treated correctly
fuck<-distance.thin(vcfR = SNPfiltR:::vcfR.example, min.distance = 10000)
fuck

#filter biallelic
#filter.biallelic()

#min mac
#this function needs work, needs to be split into min.mac function,
#and a separate function that performs repeated clustering at different missing data cutoffs
min_mac(vcfR = SNPfiltR:::vcfR.example)

#
#assess.clustering()





