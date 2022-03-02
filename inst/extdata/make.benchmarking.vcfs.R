
#read in large vcf
vcfR <- read.vcfR("~/Desktop/todiramphus.2021/populations.snps.vcf")
#subset to 100 samples and 500K SNPs
vc<-vcfR[c(1:500000),c(1:101)]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.500K.vcf.gz")

#subset to 400K
vc.400<-vc[c(1:400000),c(1:101)]
vc.400
write.vcf(x = vc.400, file="~/Desktop/benchmarking.vcfs/benchmark.400K.vcf.gz")

#300K
vc<-vc[c(1:300000),c(1:101)]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.300K.vcf.gz")

#200K
vc<-vc[c(1:200000),c(1:101)]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.200K.vcf.gz")

#100K
vc<-vc[c(1:100000),c(1:101)]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.100K.vcf.gz")

#50K
vc<-vc[c(1:50000),c(1:101)]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.50K.vcf.gz")

#20K
vc<-vc[c(1:20000),c(1:101)]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.20K.vcf.gz")

#10K
vc<-vc[c(1:10000),c(1:101)]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.10K.vcf.gz")


