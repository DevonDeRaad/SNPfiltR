
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SNPfiltR

<!-- badges: start -->

![CRAN-status](https://www.r-pkg.org/badges/version-last-release/SNPfiltR)
![CRAN-Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SNPfiltR)
![License](https://img.shields.io/badge/license-MIT-red.svg)
<!-- badges: end -->

## Installation

``` r
#Install current release from CRAN
install.packages("SNPfiltR")
```

## Citation

If you are using SNPfiltR in your pipeline for RAD analyses, I recommend
citing both SNPfiltR and the R package vcfR (which is used heavily
inside of SNPfiltR functions to read in and subset vcf files) e.g., “We
used the R packages *SNPfiltR* (DeRaad, 2021) and *vcfR* (Knaus and
Grunwald, 2017) to iteratively filter vcf files based on various quality
and missing data metrics.”

DeRaad, Devon A. 2021. Permanent DOI for SNPfiltR to come.

Knaus, Brian J., and Niklaus J. Grunwald. 2017. VCFR: a package to
manipulate and visualize variant call format data in R. Molecular
Ecology Resources 17(1):44-53.
<http://dx.doi.org/10.1111/1755-0998.12549>.

## Full documentation

To see full documentation of all functions and three detailed vignettes
illustrating use cases, please check out the [pkgdown
site](https://devonderaad.github.io/SNPfiltR/index.html) for SNPfiltR.
For a quick start, simply follow the directions below:

# Optional Step 0:

### Do quality control per sample before performing SNP calling. I have written an [RMarkdown script](https://github.com/DevonDeRaad/RADstackshelpR/blob/master/inst/extdata/fastqcr.Rmd) that uses the R package [fastqcr](https://github.com/kassambara/fastqcr) to generate a report visualizing the quality and quantity of sequencing for each sample, and recommending a subset of samples to be immediately dropped before parameter optimization (specifically useful for RADseq data). The only modification necessary for this script is the path to the folder containing the input .fastq.gz files and the path to your desired output folder. An example report generated using this script can be seen [here](https://devonderaad.github.io/RADstackshelpR/articles/quality.control.vignette.html). Because the fastq.gz files for your experiment may be large and handled remotely, an example bash script for executing this RMarkdown file as a job on a high performance computing cluster is available [here](https://github.com/DevonDeRaad/RADstackshelpR/blob/master/inst/extdata/RMarkdown.qc.submit.script.sh).

### Because this dataset was not overly large, and I used a reference based assembly that doesn’t depend on each sample contributing to the de novo building of RAD loci, I chose to skip this step and do all of the per sample quality filtering in R. I started below by reading in the vcf file using the aforementioned vcfR package:

``` r
library(SNPfiltR)
library(vcfR)
#> 
#>    *****       ***   vcfR   ***       *****
#>    This is vcfR 1.12.0 
#>      browseVignettes('vcfR') # Documentation
#>      citation('vcfR') # Citation
#>    *****       *****      *****       *****
```

``` r
#read in vcf as vcfR
vcfR <- read.vcfR("~/Desktop/aph.data/populations.snps.vcf")
```

``` r
### check the metadata present in your vcf
vcfR
#> ***** Object of Class vcfR *****
#> 115 samples
#> 87 CHROMs
#> 210,336 variants
#> Object size: 685.5 Mb
#> 57.1 percent missing data
#> *****        *****         *****

#generate popmap file. Two column popmap with the same format as stacks, and the columns must be named 'id' and 'pop'
popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 3,11))
```

# Step 1: Implement quality filters that don’t involve missing data. This is because removing low data samples will alter percentage/quantile based missing data cutoffs, so we wait to implement those until after deciding on our final set of samples for downstream analysis

### Note: I like to start with a completely unfiltered vcf file, in order to implement a cohesive, fully documentable filtering pipeline in R. However, if I have a very large unfiltered vcf file that I would like to work with in Rstudio, I sometimes perform hard filtering, or implement some conservative percentage based filters (e.g., remove all SNPs above 90% missing data) via command-line calls to the program VCFtools, in order to reduce the size of the starting vcf that I need to read into Rstudio. Just be aware that once you drop low data samples, your previously enforced (per SNP or locus) missing data % will no longer be accurate. VCFtools has faster implementations of some of the same filters implemented here, and really shines with large datasets. The R-based implementations of these filters in SNPfiltR offers interactive visualization in a cohesive R-based pipeline, but be aware that reading in excessively large files to Rstudio (typically, &gt; 1 Gb in my experience) can cause R to hang and crash. Luckily, I have found that variant sites only vcf files for most reduced-representation genomic datasets are under this size, and can be handled efficiently using Rstudio on a personal computer.

### Note: Jon Puritz has an excellent filtering tutorial that is focused specifically on filtering RADseq data: datahttps://www.ddocent.com/filtering/ Multiple functions in SNPfiltR were generated in order to follow the guidelines and suggestions laid out in e can follow these guidelines for hard filtering (he suggests minimum depth=3, gq =30), and can implement suggested filters like allele balance and max depth, here in R using SNPfiltR.

### start by visualizing the distributions of depth of sequencing and genotype quality among called genotypes, then set appropriate cutoffs for both values for this dataset.

``` r
#visualize distributions
hard_filter(vcfR=vcfR)
#> [1] "no depth cutoff provided, exploratory visualization will be generated."
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-5-2.png" width="100%" />

    #> [1] "no genotype quality cutoff provided, exploratory visualization will be generated."
    #> ***** Object of Class vcfR *****
    #> 115 samples
    #> 87 CHROMs
    #> 210,336 variants
    #> Object size: 685.5 Mb
    #> 57.1 percent missing data
    #> *****        *****         *****

    #hard filter to minimum depth of 5, and minimum genotype quality of 30
    vcfR<-hard_filter(vcfR=vcfR, depth = 5, gq = 30)
    #> [1] "32.92% of genotypes fall below a read depth of 5 and were converted to NA"
    #> [1] "2.01% of genotypes fall below a genotype quality of 30 and were converted to NA"

<img src="man/figures/README-unnamed-chunk-5-3.png" width="100%" />

Use this function to filter for allele balance From the Ddocent SNP
filtering tutorial “Allele balance: a number between 0 and 1
representing the ratio of reads showing the reference allele to all
reads, considering only reads from individuals called as heterozygous,
we expect that the allele balance in our data (for real loci) should be
close to 0.5”

### the SNPfiltR allele balance function will convert heterozygous genotypes to missing if they fall outside of the .25-.75 range.

``` r
#execute allele balance filter
vcfR<-filter_allele_balance(vcfR)
#> [1] "7.56% of het genotypes (0.39% of all genotypes) fall outside of .25 - .75 allele balance and were converted to NA"
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

Now we can execute a max depth filter (super high depth loci are likely
multiple loci stuck together into a single paralogous locus). \#\#\#
Note: we want to apply this ‘per SNP’ rather than ‘per genotype’
otherwise we will simply be removing most of the genotypes from our
deepest sequenced samples (because sequencing depth is so variable
between samples)

``` r
#visualize and pick appropriate max depth cutoff
max_depth(vcfR)
#> [1] "cutoff is not specified, exploratory visualization will be generated."
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

    #> [1] "dashed line indicates a mean depth across all SNPs of 46.7"
    #filter vcf by the max depth cutoff you chose
    vcfR<-max_depth(vcfR, maxdepth = 100)
    #> [1] "maxdepth cutoff is specified, filtered vcfR object will be returned"
    #> [1] "12.85% of SNPs were above a mean depth of 100 and were removed from the vcf"

<img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" />

Note: It may be a good idea to additionally filter out SNPs that are
significantly out of HWE if you have a really good idea of what the
population structure in your sample looks like and good sample sizes in
each pop. For this dataset, which is highly structured (many isolated
island pops) with species boundaries that are in flux, I am not going to
use a HWE filter, because I don’t feel comfortable confidently
identifying populations in which we can expect HWE. Many other programs
(such as VCFtools) can filter according to HWE if desired.

``` r
#check vcfR to see how many SNPs we have left
vcfR
#> ***** Object of Class vcfR *****
#> 115 samples
#> 78 CHROMs
#> 183,293 variants
#> Object size: 408 Mb
#> 79.79 percent missing data
#> *****        *****         *****
```

Determining which samples and SNPs to retain is alway project specific,
and is contingent on sampling, biology of the focal taxa, sequencing
idiosyncracies, etc. SNPfiltR contains functions designed to simply and
easily generate exploratory visualizations that will allow you to make
informed decisions about which samples and SNPs are of sufficient
quality to retain for downstream analyses. This is an inherently
iterative process. There is no right answer, so the best thing to do is
to look at your data, look at what some reasonable cutoffs would do to
your data, and pick one that works for you. Then as you continue to
analyze your data, make sure that filtering decisions are not driving
the patterns you see, and iteratively update your filtering approach if
you are concerned that your previous filtering choices are affecting
downstream results.

# Step 2: visualize missing data by sample. Check out the visualizations and make decision on which samples to keep for downstream analysis.

``` r
#run function to visualize samples
missing_by_sample(vcfR=vcfR, popmap = popmap)
#> Bin width defaults to 1/30 of the range of the data. Pick better value with `binwidth`.
#> Bin width defaults to 1/30 of the range of the data. Pick better value with `binwidth`.
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-9-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-9-3.png" width="100%" />

``` r
#run function to drop samples above the threshold we want from the vcf
vcfR<-missing_by_sample(vcfR=vcfR, cutoff = .91)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

    #> [1] "20 samples are above a 0.91 missing data cutoff, and were removed from VCF"
    #subset popmap to only include retained individuals
    popmap<-popmap[popmap$id %in% colnames(vcfR@gt),]

    #remove invariant sites generated by dropping individuals
    vcfR<-min_mac(vcfR, min.mac = 1)

<img src="man/figures/README-unnamed-chunk-10-2.png" width="100%" />

    #> [1] "56.14% of SNPs fell below a minor allele count of 1 and were removed from the VCF"

    #verify that missing data is not driving clustering patterns among the retained samples
    miss<-assess_missing_data_pca(vcfR=vcfR, popmap = popmap, thresholds = .8, clustering = FALSE)
    #> [1] "cutoff is specified, filtered vcfR object will be returned"
    #> [1] "76.56% of SNPs fell below a completeness cutoff of 0.8 and were removed from the VCF"
    #> Loading required namespace: adegenet
    #> Registered S3 method overwritten by 'spdep':
    #>   method   from
    #>   plot.mst ape

<img src="man/figures/README-unnamed-chunk-10-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-10-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-10-5.png" width="100%" />

``` r
#if there are still problematic samples, drop them using the following syntax
#vcfR <- vcfR[,colnames(vcfR@gt) != "A_woodhouseii_24711" & colnames(vcfR@gt) != "A_californica_45901"]
```

\#Step 3: Set the arbitrary missing data cutoff We can visualize the
effect that typical missing data cutoffs will have on both the number of
SNPs retained and the total missing data in our entire dataset. We want
to choose a cutoff that minimizes the overall missing data in the
dataset, while maximizing the total number of loci retained. \#\#\#
Note: This filter will interact with the above step, where we dropped
low data samples. A good rule of thumb is that samples shouldn’t be
above 50% missing data after applying this cutoff. So if we are
retaining low data samples out of necessity or project design, we may
have to set a more stringent cutoff at the expense of total SNPs
retained for downstream analyses. We can again use the
assess\_missing\_data\_pca() function to determine whether all retained
samples contain enough data at our chosen cutoff in order to be assigned
accurately to their species group.

``` r
#visualize missing data by SNP and the effect of various cutoffs on the missingness of each sample
missing_by_snp(vcfR)
#> [1] "cutoff is not specified, exploratory visualizations will be generated"
#> Picking joint bandwidth of 0.0645
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-11-2.png" width="100%" />

    #>    filt missingness snps.retained
    #> 1  0.30  0.30436210         47570
    #> 2  0.50  0.20138396         35398
    #> 3  0.60  0.15772868         30174
    #> 4  0.65  0.13437846         27298
    #> 5  0.70  0.11387863         24657
    #> 6  0.75  0.09326028         21836
    #> 7  0.80  0.07731504         19491
    #> 8  0.85  0.05743029         16307
    #> 9  0.90  0.03871199         12790
    #> 10 0.95  0.01867360          8013
    #> 11 1.00  0.00000000          1779

    #verify that missing data is not driving clustering patterns among the retained samples at some reasonable thresholds
    miss<-assess_missing_data_pca(vcfR=vcfR, popmap = popmap, thresholds = c(.7,.8,.85), clustering = FALSE)
    #> [1] "cutoff is specified, filtered vcfR object will be returned"

<img src="man/figures/README-unnamed-chunk-11-3.png" width="100%" />

    #> [1] "69.33% of SNPs fell below a completeness cutoff of 0.7 and were removed from the VCF"

<img src="man/figures/README-unnamed-chunk-11-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-11-5.png" width="100%" />

    #> [1] "cutoff is specified, filtered vcfR object will be returned"

<img src="man/figures/README-unnamed-chunk-11-6.png" width="100%" />

    #> [1] "76.56% of SNPs fell below a completeness cutoff of 0.8 and were removed from the VCF"

<img src="man/figures/README-unnamed-chunk-11-7.png" width="100%" /><img src="man/figures/README-unnamed-chunk-11-8.png" width="100%" />

    #> [1] "cutoff is specified, filtered vcfR object will be returned"

<img src="man/figures/README-unnamed-chunk-11-9.png" width="100%" />

    #> [1] "79.72% of SNPs fell below a completeness cutoff of 0.85 and were removed from the VCF"

<img src="man/figures/README-unnamed-chunk-11-10.png" width="100%" /><img src="man/figures/README-unnamed-chunk-11-11.png" width="100%" />

``` r
#choose a value that retains an acceptable amount of missing data in each sample, and maximizes SNPs retained while minimizing overall missing data, and filter vcf
vcfR<-missing_by_snp(vcfR, cutoff = .85)
#> [1] "cutoff is specified, filtered vcfR object will be returned"
```

<img src="man/figures/README-unnamed-chunk-11-12.png" width="100%" />

    #> [1] "79.72% of SNPs fell below a completeness cutoff of 0.85 and were removed from the VCF"

    #check how many SNPs and samples are left
    vcfR
    #> ***** Object of Class vcfR *****
    #> 95 samples
    #> 36 CHROMs
    #> 16,307 variants
    #> Object size: 111.2 Mb
    #> 5.743 percent missing data
    #> *****        *****         *****

\#Step 4: investigate the effect of a minor allele count (MAC) cutoff on
downstream inferences. MAC/MAF cutoffs can be helpful in removing
spurious and uninformative loci from the dataset, but also have the
potential to bias downstream inferences. Linck and Battey (2019) have an
excellent paper on just this topic. From the paper-

“We recommend researchers using model‐based programs to describe
population structure observe the following best practices: (a) duplicate
analyses with nonparametric methods suchas PCA and DAPC with cross
validation (b) exclude singletons (c) compare alignments with multiple
assembly parameters When seeking to exclude only singletons in
alignments with missing data (a ubiquitous problem for
reduced‐representation library preparation methods), it is preferable to
filter by the count (rather than frequency) of the minor allele, because
variation in the amount of missing data across an alignment will cause a
static frequency cutoff to remove different SFS classes at different
sites”"

Our package contains a convenient wrapper functions that streamline
investigation of the effects of various MAC cutoffs.

``` r
#investigate clustering patterns with and without a minor allele cutoff
#use min.mac() to investigate the effect of multiple cutoffs
vcfR.mac<-min_mac(vcfR = vcfR, min.mac = 2)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

    #> [1] "30.87% of SNPs fell below a minor allele count of 2 and were removed from the VCF"

    #assess clustering without MAC cutoff
    miss<-assess_missing_data_tsne(vcfR, popmap, clustering = FALSE)

<img src="man/figures/README-unnamed-chunk-12-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-12-3.png" width="100%" />

``` r
#assess clustering with MAC cutoff
miss<-assess_missing_data_tsne(vcfR.mac, popmap, clustering = FALSE)
```

<img src="man/figures/README-unnamed-chunk-12-4.png" width="100%" /><img src="man/figures/README-unnamed-chunk-12-5.png" width="100%" />

``` r

#based on these visualizations, singletons are not obviously biasing clustering patterns, so I will leave them in for now. If I want to run a program like STRUCTURE, where singletons are known to bias inference, I can write out the vcf with singletons removed as well:
#vcfR::write.vcf(vcfR.mac, file = "~/Downloads/scrub.jay.mac.vcf.gz")
```

### Note: the function vcfR::write.vcf() automatically writes a gzipped vcf file, so be sure to add the suffix .gz to the name of your output file.

### Make sure that the depth and genotype quality look consistent across SNPs and samples

``` r
#plot depth per snp and per sample
dp <- extract.gt(vcfR, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

``` r
#plot genotype quality per snp and per sample
gq <- extract.gt(vcfR, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)
```

<img src="man/figures/README-unnamed-chunk-13-2.png" width="100%" />

### Write out the filtered vcf, and a linkage filtered version, for downstream analyses

``` r
#fix mislabeled sample
colnames(vcfR@gt)[colnames(vcfR@gt) == "A_californica_334171"]<-"A_woodhouseii_334171"
#write out vcf with all SNPs
vcfR::write.vcf(vcfR, "~/Downloads/aphelocoma.filtered.vcf.gz")

#linkage filter vcf to thin SNPs to one per 500bp
vcfR.thin<-distance_thin(vcfR, min.distance = 500)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  81%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
#> [1] "2803 out of 16307 input SNPs were not located within 500 base-pairs of another SNP and were retained despite filtering"

#write out thinned vcf
vcfR::write.vcf(vcfR.thin, "~/Downloads/aphelocoma.filtered.thinned.vcf.gz")
```
