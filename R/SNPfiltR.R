#' SNPfiltR: A package for interactively visualizing and filtering SNP datasets
#'
#' The SNPfiltR package allows users to interactively visualize the effects of
#' relevant filters on their datasets in order to optimize filtering parameters
#' rather than simply following historical precedent. Each function takes a
#' variant call format (vcf) file, stored in local memory as a vcfR object, as input.
#' Most functions can be run without specified cutoffs, in order to visualize the
#' distribution of the parameter of interest in your given dataset. Then the same
#' function can be run with a specified cutoff, and a filtered vcfR object
#' will be returned. For detailed documentation and vignettes showing fully implemented
#' SNP filtering pipelines, please go to: devonderaad.github.io/SNPfiltR
#'
#' @docType package
#' @name SNPfiltR
NULL
#> NULL
