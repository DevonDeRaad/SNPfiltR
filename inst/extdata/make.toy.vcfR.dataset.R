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
vcfR.example<-vcffog[sample.int(151015, 500),]

#save this vcf
write.vcf(vcfR.example, file = "~/Downloads/SNPfilter.example.vcf.gz")

#make example popmap for example vcfR
popmap<-data.frame(id=colnames(SNPfiltR:::vcfR.example@gt)[-1], pop=sub(".*H_", "", colnames(SNPfiltR:::vcfR.example@gt)[-1]))

#save this popmap as available data for the SNPfiltR package
#write.csv(popmap, "")

#save this vcfR as available data for the SNPfiltR package
usethis::use_data(vcfR.example, internal = TRUE, overwrite = T)

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
filter_biallelic(vcfR = SNPfiltR:::vcfR.example)
vc<-read.vcfR("/Users/devder/Dropbox/Emily_Towhees/towhees/towhees_genomic/Pipmac2020.filtered.vcf")
f<-filter_biallelic(vcfR = vc)
filter_biallelic(vcfR = f)


#min mac
#this function now only does mac filtering
min_mac(vcfR = SNPfiltR:::vcfR.example)
f<-min_mac(vcfR = SNPfiltR:::vcfR.example, min.mac = 1)
min_mac(vcfR = f, min.mac = 4)

#assess clustering
assess_clustering(vcfR = SNPfiltR:::vcfR.example,
                  popmap = popmap,
                  thresholds = c(.6,.7,.8,.9))

#
assess_missing_data(vcfR = SNPfiltR:::vcfR.example,
                  popmap = popmap,
                  thresholds = c(.6,.7,.8,.9))




