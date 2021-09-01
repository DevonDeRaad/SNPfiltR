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
#usethis::use_data(vcfR.example, internal = TRUE, overwrite = T)

#clean global environment
rm(list=ls())

#load and test examples

#missing by snp
missing_by_snp(vcfR = SNPfiltR:::vcfR.example)
f<-missing_by_snp(vcfR = SNPfiltR:::vcfR.example, cutoff = .6)
missing_by_snp(vcfR = f, cutoff = .6)

#missing by sample
missing_by_sample(vcfR = SNPfiltR:::vcfR.example)
f<-missing_by_sample(vcfR = SNPfiltR:::vcfR.example, cutoff = .7)
missing_by_sample(vcfR = f, cutoff = .7)

#hard filter
f<-hard_filter(vcfR = SNPfiltR:::vcfR.example, depth = 5, gq = 30)
hard_filter(vcfR = f, depth = 5, gq = 30)

#max depth
max_depth(vcfR = SNPfiltR:::vcfR.example)
f<-max_depth(vcfR = SNPfiltR:::vcfR.example, maxdepth = 100)
max_depth(vcfR = f, maxdepth = 100)

#allele balance
#fix histograms
f<-filter_allele_balance(vcfR = SNPfiltR:::vcfR.example)
filter_allele_balance(vcfR = f)

#distance thin
#function now works correctly
f<-distance_thin(vcfR = SNPfiltR:::vcfR.example, min.distance = 1000)
distance_thin(vcfR = f, min.distance = 1000)

#filter biallelic
#filter_biallelic()

#min mac
#this function now only does mac filtering
min_mac(vcfR = SNPfiltR:::vcfR.example)
f<-min_mac(vcfR = SNPfiltR:::vcfR.example, min.mac = 4)
min_mac(vcfR = f, min.mac = 4)

#
#assess_clustering()





