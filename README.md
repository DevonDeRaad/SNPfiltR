
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

## Overview

This is a basic example which shows you how to solve a common problem:
