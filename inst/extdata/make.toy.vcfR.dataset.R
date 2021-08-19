#devtools::install_github("DevonDeRaad/SNPfiltR")

#clean package
#library(pkgdown)
#clean_site(pkg = "SNPfiltR")
#library(devtools)
#clean_vignettes(pkg = "SNPfiltR")

#read in vcf and make subset example datasets
library(vcfR)

vcffog<-read.vcfR("~/Desktop/hipposideros/n3.vcf")

vc.500<-vcffog[sample.int(151015, 500), c(1:20)]

write.vcf(vc.500, file = "~/Desktop/SNPfiltR/inst/extdata/unfiltered.vcf.gz")

vcfR<-read.vcfR("~/Desktop/SNPfiltR/inst/extdata/unfiltered.vcf.gz")
usethis::use_data(vcfR)

missing.by.sample(vcfR = system.file("extdata", "unfiltered.vcf.gz", package = "RADstackshelpR"))



#Now do big M
vc.1<-vcffog[sample.int(110576, 100), c(1:20)]
vc.2<-vcffog[sample.int(110576, 90), c(1:20)]
vc.3<-vcffog[sample.int(110576, 80), c(1:20)]
vc.4<-vcffog[sample.int(110576, 70), c(1:20)]
vc.5<-vcffog[sample.int(110576, 60), c(1:20)]
vc.6<-vcffog[sample.int(110576, 50), c(1:20)]
vc.7<-vcffog[sample.int(110576, 40), c(1:20)]
vc.8<-vcffog[sample.int(110576, 30), c(1:20)]

#write.vcf(vc.1, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM1.vcf.gz")
#write.vcf(vc.2, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM2.vcf.gz")
#write.vcf(vc.3, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM3.vcf.gz")
#write.vcf(vc.4, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM4.vcf.gz")
#write.vcf(vc.5, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM5.vcf.gz")
#write.vcf(vc.6, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM6.vcf.gz")
#write.vcf(vc.7, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM7.vcf.gz")
#write.vcf(vc.8, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM8.vcf.gz")

opt.bigm<- optimize_bigM(M1 = system.file("extdata", "bigM1.vcf.gz", package = "RADstackshelpR"),
                         M2 = system.file("extdata", "bigM2.vcf.gz", package = "RADstackshelpR"),
                         M3 = system.file("extdata", "bigM3.vcf.gz", package = "RADstackshelpR"),
                         M4 = system.file("extdata", "bigM4.vcf.gz", package = "RADstackshelpR"),
                         M5 = system.file("extdata", "bigM5.vcf.gz", package = "RADstackshelpR"),
                         M6 = system.file("extdata", "bigM6.vcf.gz", package = "RADstackshelpR"),
                         M7 = system.file("extdata", "bigM7.vcf.gz", package = "RADstackshelpR"),
                         M8 = system.file("extdata", "bigM8.vcf.gz", package = "RADstackshelpR"))

vis_loci(opt.bigm)

#now do n
#write.vcf(vc.8, file = "~/Desktop/RADstackshelpR/inst/extdata/nequalsmminus1.vcf.gz")
#write.vcf(vc.6, file = "~/Desktop/RADstackshelpR/inst/extdata/nequalsm.vcf.gz")
#write.vcf(vc.7, file = "~/Desktop/RADstackshelpR/inst/extdata/nequalsmplus1.vcf.gz")

opt.n<- optimize_n(nequalsMminus1 = system.file("extdata", "nequalsmminus1.vcf.gz", package = "RADstackshelpR"),
                   nequalsM = system.file("extdata", "nequalsm.vcf.gz", package = "RADstackshelpR"),
                   nequalsMplus1 = system.file("extdata", "nequalsmplus1.vcf.gz", package = "RADstackshelpR"))

vis_loci(opt.n)


