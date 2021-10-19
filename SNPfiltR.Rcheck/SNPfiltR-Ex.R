pkgname <- "SNPfiltR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SNPfiltR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("assess_missing_data_pca")
### * assess_missing_data_pca

flush(stderr()); flush(stdout())

### Name: assess_missing_data_pca
### Title: Vizualise how missing data thresholds affect sample clustering
### Aliases: assess_missing_data_pca

### ** Examples

assess_missing_data_pca(vcfR = SNPfiltR::vcfR.example,
popmap = SNPfiltR::popmap,
thresholds = c(.6,.8))



cleanEx()
nameEx("assess_missing_data_tsne")
### * assess_missing_data_tsne

flush(stderr()); flush(stdout())

### Name: assess_missing_data_tsne
### Title: Vizualise how missing data thresholds affect sample clustering
### Aliases: assess_missing_data_tsne

### ** Examples

assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
popmap = SNPfiltR::popmap,
thresholds = .8)



cleanEx()
nameEx("distance_thin")
### * distance_thin

flush(stderr()); flush(stdout())

### Name: distance_thin
### Title: Filter a vcf file based on distance between SNPs on a given
###   scaffold
### Aliases: distance_thin

### ** Examples

distance_thin(vcfR = SNPfiltR::vcfR.example, min.distance = 1000)



cleanEx()
nameEx("filter_allele_balance")
### * filter_allele_balance

flush(stderr()); flush(stdout())

### Name: filter_allele_balance
### Title: Filter out heterozygous genotypes failing an allele balance
###   check
### Aliases: filter_allele_balance

### ** Examples

filter_allele_balance(vcfR = SNPfiltR::vcfR.example)



cleanEx()
nameEx("filter_biallelic")
### * filter_biallelic

flush(stderr()); flush(stdout())

### Name: filter_biallelic
### Title: Remove SNPs with more than two alleles
### Aliases: filter_biallelic

### ** Examples

filter_biallelic(vcfR = SNPfiltR::vcfR.example)



cleanEx()
nameEx("hard_filter")
### * hard_filter

flush(stderr()); flush(stdout())

### Name: hard_filter
### Title: Hard filter a vcf file by depth and genotype quality (gq)
### Aliases: hard_filter

### ** Examples

hard_filter(vcfR = SNPfiltR::vcfR.example, depth = 5)
hard_filter(vcfR = SNPfiltR::vcfR.example, depth = 5, gq = 30)



cleanEx()
nameEx("max_depth")
### * max_depth

flush(stderr()); flush(stdout())

### Name: max_depth
### Title: Vizualise and filter based on mean depth across all called SNPs
### Aliases: max_depth

### ** Examples

max_depth(vcfR = SNPfiltR::vcfR.example)
max_depth(vcfR = SNPfiltR::vcfR.example, maxdepth = 100)



cleanEx()
nameEx("min_mac")
### * min_mac

flush(stderr()); flush(stdout())

### Name: min_mac
### Title: Vizualise, filter based on Minor Allele Count (MAC)
### Aliases: min_mac

### ** Examples

min_mac(vcfR=SNPfiltR::vcfR.example)



cleanEx()
nameEx("missing_by_sample")
### * missing_by_sample

flush(stderr()); flush(stdout())

### Name: missing_by_sample
### Title: Vizualise missing data per sample, remove samples above a
###   missing data cutoff
### Aliases: missing_by_sample

### ** Examples

missing_by_sample(vcfR = SNPfiltR::vcfR.example)
missing_by_sample(vcfR = SNPfiltR::vcfR.example, cutoff = .7)



cleanEx()
nameEx("missing_by_snp")
### * missing_by_snp

flush(stderr()); flush(stdout())

### Name: missing_by_snp
### Title: Vizualise missing data per SNP, remove SNPs above a missing data
###   cutoff
### Aliases: missing_by_snp

### ** Examples

missing_by_snp(vcfR = SNPfiltR::vcfR.example)
missing_by_snp(vcfR = SNPfiltR::vcfR.example, cutoff = .6)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
