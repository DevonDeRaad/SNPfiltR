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

`SNPfiltR` is an R [@R Core Team:2019] package for interactively and reproducibly filtering SNP (single nucleotide polymorphism) data.

Restriction-site Associated DNA sequencing (RADseq) is a commonly used method for rapidly and cost-effectively generating sequence data from thousands of genome-wide loci. One great advantage of this approach is that these genome-wide loci can be assembled de novo (i.e., without aligning the sequences to a reference genome) making RAD sequencing ideal for non-model systems for which a closely-related reference genome is not available. Additionally, RADseq offers reasonable file-sizes and cost of sequencing, which has resulted in the rapid adoption of this convenient approach for answering questions in both population genetics and phylogenetics over the last decade.

Unfortunately, RADseq data also comes with the downside of a dizzying array of computational choices which must be made during bioinformatic processing, especially during de novo locus assembly and single nucleotide polymorphism (SNP) calling. Rampant missing data, uneven sequencing amongst samples, idiosyncratic features of the focal genome, and stochastic factors associated with a given sequencing run, can all contribute to making the optimal parameters for assembling loci and calling SNPs from a given RAD sequencing run a moving target. The software pipeline Stacks [@Catchen:2011] is designed specifically for assembling loci and calling SNPs from RADseq data, and is designed modularly, in order to give end-users control over a variety of assembly parameters throughout the pipeline. The paper 'Lost in Parameter Space' [@Paris:2018] performed detailed explorations of parameter space for the de novo assembly of RAD loci from three empirical datasets, resulting in a clearly defined set of best practices for de novo assembly of RAD loci using the Stacks pipeline.

# Statement of need

A number of R packages are dedicated to manipulating SNP data, but most focus either on reading, writing, and directly manipulating SNP data (e.g., vcfR [@Knaus:2017] and dartR [cite]), or performing downstream phylogenetic and population genetic analyses (e.g., APE [cite], pegas [cite], stAMPP [cite], SNPrelate [cite], adegenet [cite], and introgress [cite]). In contrast, `SNPfiltR` is specifically designed for filtering SNP data, especially for reduced-representation genomic datasets, by providing wrapper scripts to perform exploratory data visualizations and directly filter SNPs and genotypes based on metrics such as quality, depth of coverage, and amount of missing data. `SNPfiltR` makes the SNP filtering process interactive, reproducible, and easily documentable thanks to the robust suite of open access resources available for computing in R (specifically Rstudio [cite] for interactivity and Rmarkdown [cite] for reproducibility). Users will note that many of the filters implemented by `SNPfiltR` are already available elsewhere, specifically via command line programs like VCFtools [cite], and GATK [cite]. We show that `SNPfiltR` is comparable in performance with these command-line programs when run on small to moderately sized SNP datasets, while larger SNP datasets are hindered by the memory requirements associated with reading large files into R (Fig. _)

`RADstackshelpR` contains a series of wrapper functions which make use internally of the R package vcfR [@Knaus:2017] to read in variant call format (vcf) files which are output directly by Stacks, and return a data frame detailing the number of polymorphic loci and SNPs contained in each vcf file at an 80% completeness cutoff ('R80') [@Paris:2018]. These output data frames can then be input directly into convenient visualization functions which internally utilize the R packages ggplot2 [@Wickham:2020], ggridges [@Claus:2018], and gridExtra [@Auguie:2017] to quickly and reproducibly generate publication quality figures detailing the optimization process for your dataset.

Now, rather than relying on series' of spreadsheets and disparate script files to haphazardly explore different parameter combinations, biologists looking to tackle de novo assembly of RADseq data can follow an explicit pipeline for thoroughly and repeatably exploring parameter space according to best practices and visualizing their results for publication.

![Fig 1. Multi-panel figure visualizing the entire optimization process facilitated by `RADstackshelpR`.\label{fig:overview}](fig1.png)

# Acknowledgements
I would like to acknowledge the authors of the R package vcfR, Brian J. Knaus, and Nikolas J. Grünwald, for developing and documenting resources for reading and manipulating SNP data using R, specifically the exceptional GitHub site: https://knausb.github.io/vcfR_documentation/index.html which provided invaluable guidance in the development of this package.

# References
