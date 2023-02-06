#' Remove SNPs with more than two alleles
#'
#' This function simply removes any SNPs from the vcf file which contains more than
#' two alleles. Many downstream applications require
#' SNPs to be biallelic, so this filter is generally a good idea during processing.
#'
#' @param vcfR a vcfR object
#' @return a vcfR object with SNPs containing more than two alleles removed
#' @examples
#' filter_biallelic(vcfR = SNPfiltR::vcfR.example)
#' @export
filter_biallelic <- function(vcfR){

    #if vcfR is not class vcfR, fail gracefully
    if (class(vcfR) != "vcfR"){
      stop("specified vcfR object must be of class 'vcfR'")
    }

    #if "ALT" is not specified, throw an informative error
    if ("ALT" %in% colnames(vcfR@fix) == FALSE){
      #if all "ALT" values are null, exit
      stop(" 'ALT' column not specified in the input vcf, can't determine number of alternate alleles")
    }

    #if "ALT" is not informative, throw an informative error
    if (sum(is.na(vcfR@fix[,"ALT"])) == length(vcfR@fix[,"ALT"])){
      #if all "ALT" values are null, exit
      stop("informative values are not supplied in the 'ALT' column of the input vcf, can't determine number of alternate alleles")
    }

    #store vector of number of alternate alleles at each SNP
    v<-nchar(gsub(",","",vcfR@fix[,"ALT"]))
    #add 1 to each value to account for the reference allele
    v<-v+1

    #make histogram showing the distribution of number of alleles across all SNPs
    graphics::hist(v,
         main= "distribution of alleles present in vcf",
         xlab= "number of alleles",
         ylab= "number of SNPs",
         xlim= c(.5,4.5))
    graphics::abline(v=2.5,
         col="red")

    #number of input SNPs
    p<-nrow(vcfR@fix)

    #number of SNPs to be removed (based on ALT column)
    q<-p-sum(nchar(vcfR@fix[,"ALT"]) == 1)

    #% of SNPs that will be removed
    y<-q/p

    #print for user
    message(q, " SNPs, ", round(y, 3), "% of all input SNPs, contained more than 2 alleles, and were removed from the VCF")

    #filter the vcfR based on the "ALT" column
    vcfR<-vcfR[nchar(vcfR@fix[,"ALT"]) == 1,]

    #return vcf
    return(vcfR)

}

