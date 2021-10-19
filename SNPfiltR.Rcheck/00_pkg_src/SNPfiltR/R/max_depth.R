#' Vizualise and filter based on mean depth across all called SNPs
#'
#' This function can be run in two ways: 1) specify vcfR object only. This will
#' visualize the distribution of mean depth per sample across all SNPs in your vcf file,
#' and will not alter your vcf file. 2) specify vcfR object, and set 'maxdepth' = 'integer value'.
#' This option will show you where your specified cutoff falls in the distribution of SNP depth,
#' and remove all SNPs with a mean depth above the specified threshold from the vcf.
#' Super high depth loci are likely multiple loci stuck together into a single paralogous locus.
#' Note: This function filters on a 'per SNP' basis rather than a 'per genotype' basis,
#' otherwise it would disproportionately remove genotypes from our deepest sequenced samples
#' (because sequencing depth is so variable between samples).
#'
#' @param vcfR a vcfR object
#' @param maxdepth an integer specifying the maximum mean depth for a SNP to be retained
#' @return The vcfR object input, with SNPs above the 'maxdepth' cutoff removed
#' @examples
#' max_depth(vcfR = SNPfiltR::vcfR.example)
#' max_depth(vcfR = SNPfiltR::vcfR.example, maxdepth = 100)
#' @export
max_depth <- function(vcfR,
                      maxdepth=NULL){

  #if specified vcfR is not class 'vcfR', fail gracefully
  if (class(vcfR) != "vcfR"){
    stop("specified vcfR object must be of class 'vcfR'")
  }

  #if maxdepth cutoff is not specified, start here
  if (is.null(maxdepth)){

    #print user message
    print("cutoff is not specified, exploratory visualization will be generated.")

    #extract depth from the vcf
    dp.matrix<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)

    #calculate vector of depth per SNP
    snpdepth<-rowSums(dp.matrix, na.rm = TRUE)/rowSums(is.na(dp.matrix) == FALSE)

    #set plotting parameters
    graphics::par(mfrow=c(2,1))
    #plot histogram of depth
    graphics::hist(snpdepth,
         xlab = "mean of the depth of all samples successfully genotyped at a given SNP",
         main="histogram showing the depth of all called SNPs")
    graphics::abline(v=mean(snpdepth, na.rm = TRUE),
           col="red",
           lty="dashed")

    #zoomed in histogram
    graphics::hist(snpdepth[snpdepth < 200],
         xlab = "mean of the depth of all samples successfully genotyped at a given SNP",
         main ="distribution of SNPs below a depth of 200")
    graphics::abline(v=mean(snpdepth, na.rm = TRUE),
           col="red",
           lty="dashed")

    #print
    print(paste0("dashed line indicates a mean depth across all SNPs of ",round(mean(snpdepth, na.rm = TRUE),1)))

  }

  #if the maxdepth cutoff has been specified, start here
  else {

    if (is.numeric(maxdepth) != "TRUE"){
      stop("specified max depth cutoff must be numeric")
    }

    #print user message
    print("maxdepth cutoff is specified, filtered vcfR object will be returned")

    #extract depth from the vcf
    dp.matrix<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)

    #calculate vector of depth per SNP
    snpdepth<-rowSums(dp.matrix, na.rm = TRUE)/rowSums(is.na(dp.matrix) == FALSE)

    #plot the maxdepth cutoff
    graphics::hist(snpdepth,
         xlab = "mean of the depth of all samples successfully genotyped at a given SNP",
         main ="max depth cutoff")
    graphics::abline(v=maxdepth,
           col="red")

    #calculate % of snps that fail the max depth filter
    i<-round(sum(snpdepth > maxdepth, na.rm = TRUE)/length(snpdepth)*100, 2)

    print(paste0(i, "% of SNPs were above a mean depth of ", maxdepth, " and were removed from the vcf"))

    #filter vcf
    vcfR<-vcfR[snpdepth < maxdepth,]

    #return vcfR
    return(vcfR)
  #close else statement
  }
#close function
}

