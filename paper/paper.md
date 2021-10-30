---
title: 'SNPfiltR: an R package for interactive and reproducible SNP filtering'
tags:
  - R
  - single nucleotide polymorphisms
  - SNP filtering
  - bioinformatics
  - reproducibility
authors:
  - name: Devon A. DeRaad^[Corresponding author]
    orcid: 0000-0003-3105-985X
    affiliation: 1
affiliations:
 - name: Department of Ecology & Evolutionary Biology and Biodiversity Institute, University of Kansas, Lawrence, Kansas, USA
   index: 1
date: 22 October 2021
bibliography: paper.bib
---

# Summary

`SNPfiltR` is an R [@R Core Team:2019] package that facilitates the development of customizable, interactive, and reproducible SNP (single nucleotide polymorphism) filtering pipelines using the language and tools provided by the open source R community.

As next-generation (i.e. massively parallel, short read) sequencing has become the ubiquitous approach for population genetic and phylogenetic investigations, SNP datasets have enjoyed a similar rise in popularity. SNP datasets typically contain thousands to millions of base-calls for each individual sample (genotypes), located at thousands to millions of variable sites (SNPs) throughout the genome of the focal taxa. Depending on the sequencing approach, these genotypes and SNPs may suffer from a variety of issues such as sequencing errors, poor sequencing quality, and missing data, all of which should be addressed before performing downstream population genetic and phylogenetic analyses. Furthermore, individual samples might suffer from low sequencing coverage, contamination, or lack of homology among sampled fragments, resulting in individual samples that cannot be confidently used in downstream analyses, and need to be dropped from the dataset.

As reproducibility is becoming more widely acknowledged as critical for the future of science, a major challenge continues to be effectively documenting the iterative and often patchwork process of investigating, visualizing, and cleaning large datasets to prepare them for final analyses that will be presented in publication format. Many pipelines exist for filtering SNP data, and most investigators employ a combination of approaches tailored to their individual sequencing method and to each specific dataset. Many powerful command-line based tools such as VCFtools [cite], GATK [cite], STACKS [cite], Ddocent [cite], and ANGSD [cite] exist for the purpose of filtering SNPs and genotypes based on various quality metrics. While these programs are highly efficient, especially for massive datasets, the command-line interface does not lend itself easily to graphical visualization, and as a result many investigators rely on series of homebrewed scripts for investigating their data, potentially including poorly documented and non-reproducible (i.e., requiring manual input) steps within their pipelines. One common approach involves generating summary statistic files using one or more of these command-line based methods, to be visualized and investigated using the powerful data visualization tools offered by the R computing language.

`SNPfiltR` is designed to bridge these functionalities, while allowing for the needs of customizability, interactivity, and reproducibility in any SNP filtering pipeline. The `SNPfiltR` package bundles a suite of convenient wrapper functions which can be used to investigate, visualize, and filter moderately sized SNP datasets, rapidly and interactively, within a single R session. The package is publicly available on CRAN via the syntax: install.packages('SNPfiltR') and is thoroughly documented including vignettes implementing complete example pipelines for filtering SNP data derived from both RADseq <https://devonderaad.github.io/SNPfiltR/articles/scrub-jay-RADseq-vignette.html> and UCE <https://devonderaad.github.io/SNPfiltR/articles/scrub-jay-UCE-vignette.html> sequencing approaches.

# Statement of need

A number of R packages are dedicated to manipulating SNP data, but most focus either on direct manipulation of SNP data (e.g., vcfR [@Knaus:2017] and dartR [cite]), or on performing downstream phylogenetic and population genetic analyses (e.g., APE [cite], pegas [cite], stAMPP [cite], SNPrelate [cite], adegenet [cite], and introgress [cite]). `SNPfiltR` is designed as a bridge between these functionalities, leveraging the ability of the vcfR package to efficiently read vcf files into an R working environment as 'vcfR' objects, to rapidly and reproducibly filter and return vcfR objects, which can be exported in standard vcf format for downstream analyses. `SNPfiltR` is designed to make filtering SNP data simple, interactive, and reproducible, by building on functions from the R package vcfR [@Knaus:2017] which can efficiently read and store vcf files into an R working environment as objects of class 'vcfR'. All `SNPfiltR` functions take vcfR objects as input, which can be efficiently visualized and filtered based on an array of quality metrics, and return vcfR objects as output, which can then be written to local directories in standard vcf format using the vcfR function vcfR::write.vcf(). Because `SNPfiltR` relies on the vcfR package for both input formatting and writing output vcfR objects in standard vcf format for downstream analyses, we suggest that users cite vcfR alongside SNPfiltR when using this package as part of a SNP filtering pipeline.

`SNPfiltR` is especially designed for reduced-representation genomic datasets, for which unfiltered vcf files containing all called SNPs are generally less than a Gigabyte in total size. These moderately sized input files are highly computationally tractable, and can be  usiby providing wrapper scripts to perform exploratory data visualizations and directly filter SNPs and genotypes based on metrics such as quality, depth of coverage, and amount of missing data. `SNPfiltR` makes the SNP filtering process interactive, reproducible, and easily documentable thanks to the robust suite of open access resources available for computing in R (specifically Rstudio [cite] for interactivity and Rmarkdown [cite] for reproducibility). Users will note that many of the filters implemented by `SNPfiltR` are already available elsewhere, specifically via command line programs like VCFtools [cite], and GATK [cite]. We show that `SNPfiltR` is comparable in performance with these command-line programs when run on small to moderately sized SNP datasets, while larger SNP datasets are hindered by the memory requirements associated with reading large files into R (Fig. _)

The `SNPfiltR` package provides to users the novel opportunity to follow an interactive and customizable SNP filtering pipeline, fully implemented in R, therefore fully reproducible and documentable. Performance testing indicates that the SNP filtering functions offered by `SNPfiltR` are comparable in performance on small to moderate sized genomic datasets, even compared with powerful, highly-efficient command-line based programs designed for SNP filtering (Fig. 1). As an R based implementation, `SNPfiltR` struggles to scale with large genomic datasets, where local memory availability becomes a limiting factor. Nonetheless, we show that for larger genomic datasets, `SNPfiltR` can be readily integrated into a reproducible pipeline alongside aforementioned powerful, highly-efficient command-line based programs designed for SNP filtering, such as VCFtools. Ultimately, the universe of elegant, open-source R based tools such as Rstudio and Rmarkdown are ideal for facilitating both interactivity and reproducibility, and ought to lead bioinformaticians to consider an R based pipeline for streamlining the often complicated and iterative process of optimizing filtering parameters for next-generation sequencing datasets. Now, `SNPfiltR` offers a suite of powerful, modularized functions for SNP visualization and filtering, allowing users to build custom, interactive, reproducible SNP filtering pipelines, all using the widely adopted R programming language.

![Fig 1. Multi-panel figure visualizing the entire optimization process facilitated by `RADstackshelpR`.\label{fig:overview}](fig1.png)

# Acknowledgements
I would like to acknowledge the authors of the R package vcfR, Brian J. Knaus, and Nikolas J. Grünwald, for developing and documenting resources for reading and manipulating SNP data using R, specifically the exceptional GitHub site: https://knausb.github.io/vcfR_documentation/index.html which provided invaluable guidance in the development of this package.

# References
