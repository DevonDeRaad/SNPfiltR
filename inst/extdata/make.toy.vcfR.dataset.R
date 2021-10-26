#devtools::install_github("DevonDeRaad/SNPfiltR")

#clean package
#library(pkgdown)
#clean_site(pkg = "SNPfiltR")
#library(devtools)
#clean_vignettes(pkg = "SNPfiltR")

#load packages
library(vcfR)
library(SNPfiltR)

#read in vcf and make subset example datasets
vcffog<-read.vcfR("~/Desktop/aph.data/populations.snps.vcf")

#subset to 20 samples (5 cali, 5 island, 5 wood, 5 florida)
vcffog<-vcffog[,c(1:6,49:53,77:81,41:45)]

#remove invariant SNPs
vcffog<-min_mac(vcffog, min.mac = 1)

#subset to 500 random SNPs
vcfR.example<-vcffog[sample.int(2725, 500),]

#save this vcf
#write.vcf(vcfR.example, file = "~/Downloads/SNPfilter.example.vcf.gz")

#make example popmap for example vcfR
popmap<-data.frame(id=colnames(vcfR.example@gt)[-1], pop=gsub("A_", "", colnames(vcfR.example@gt)[-1]))
popmap$pop<-gsub("_.*", "", popmap$pop)
popmap$pop<-as.factor(popmap$pop)

#save this popmap as available data for the SNPfiltR package
#usethis::use_data(popmap,vcfR.example, internal = FALSE, overwrite = T)

#save this vcfR as available data for the SNPfiltR package
#usethis::use_data(vcfR.example, internal = TRUE, overwrite = F)

#clean global environment
#rm(list=ls())

#load and test examples

#missing by snp
missing_by_snp(vcfR = SNPfiltR::vcfR.example)
f<-missing_by_snp(vcfR = SNPfiltR::vcfR.example, cutoff = .6)
missing_by_snp(vcfR = f, cutoff = .6)

#missing by sample
missing_by_sample(vcfR = SNPfiltR::vcfR.example)
f<-missing_by_sample(vcfR = SNPfiltR::vcfR.example, cutoff = .7)
missing_by_sample(vcfR = f, cutoff = .7)

#hard filter
hard_filter(vcfR = SNPfiltR::vcfR.example)
f<-hard_filter(vcfR = SNPfiltR::vcfR.example, depth = 5, gq = 30)
hard_filter(vcfR = f, depth = 5, gq = 30)

#max depth
max_depth(vcfR = SNPfiltR::vcfR.example)
f<-max_depth(vcfR = SNPfiltR::vcfR.example, maxdepth = 100)
max_depth(vcfR = f, maxdepth = 100)

#allele balance
#fix histograms
f<-filter_allele_balance(vcfR = SNPfiltR::vcfR.example)
filter_allele_balance(vcfR = f)

#distance thin
#function now works correctly
f<-distance_thin(vcfR = SNPfiltR::vcfR.example, min.distance = 1000)
distance_thin(vcfR = f, min.distance = 1000)

#filter biallelic
filter_biallelic(vcfR = SNPfiltR::vcfR.example)
#vc<-read.vcfR("/Users/devder/Dropbox/Emily_Towhees/towhees/towhees_genomic/Pipmac2020.filtered.vcf")
#f<-filter_biallelic(vcfR = vc)
#filter_biallelic(vcfR = f)


#min mac
#this function now only does mac filtering
min_mac(vcfR = SNPfiltR::vcfR.example, min.mac = 2)

#assess clustering
#
assess_missing_data_pca(vcfR = SNPfiltR::vcfR.example,
                  popmap = SNPfiltR::popmap,
                  thresholds = c(.6,.8))

#
assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
                    popmap = SNPfiltR::popmap,
                    thresholds = .9)



