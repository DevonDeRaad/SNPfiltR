
#read in large vcf
vcfR <- read.vcfR("~/Desktop/todiramphus.2021/populations.snps.vcf")
#subset to 100 samples and 500K SNPs
vc<-vcfR[c(1:500000),c(1:101)]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.500K.vcf.gz")

#subset to 400K
vc<-vc[sample(c(1:500000), 400000),]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.400K.vcf.gz")

#300K
vc<-vc[sample(c(1:400000), 300000),]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.300K.vcf.gz")

#200K
vc<-vc[sample(c(1:300000), 200000),]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.200K.vcf.gz")

#100K
vc<-vc[sample(c(1:200000), 100000),]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.100K.vcf.gz")

#50K
vc<-vc[sample(c(1:100000), 50000),]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.50K.vcf.gz")

#20K
vc<-vc[sample(c(1:50000), 20000),]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.20K.vcf.gz")

#10K
vc<-vc[sample(c(1:20000), 10000),]
vc
write.vcf(x = vc, file="~/Desktop/benchmarking.vcfs/benchmark.10K.vcf.gz")


