#' Popmap for example scrub-jay vcfR file
#'
#' A dataset containing the sample name and population assignment for the 20
#'  scrub-jay samples in SNPfilR::vcfR.example . The variables are as follows:
#'
#' \itemize{
#'   \item id. unique sample identifier
#'   \item pop. population assignment for each individual
#'   }
#' @docType data
#' @keywords datasets
#' @name popmap
#' @usage data(popmap)
#' @format A data frame with 20 rows and 2 variables
"popmap"

#' Example scrub-jay vcfR file
#'
#' A vcfR object containing 500 SNPs for 20 individuals. Species assignments for
#' each individual can be accessed via SNPfiltR::popmap
#' @docType data
#' @keywords datasets
#' @name vcfR.example
#' @usage data(vcfR.example)
#' @format A vcfR object containing 500 SNPs for 20 individuals
"vcfR.example"
