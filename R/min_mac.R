#' Vizualise, filter based on Minor Allele Count (MAC)
#'
#' This function can be run in two ways: 1) Without 'min.mac' specified. This will return
#' a folded site frequency spectrum (SFS), without performing any filtering on the vcf file. Or
#' 2) With 'min.mac' specified. This will also print the folded SFS and
#' show you where your specified min. MAC count falls. It will then return your vcfR object
#' with SNPs falling below your min. MAC threshold removed.
#' Note: previous filtering steps (especially removing samples) may have resulted
#' in invariant SNPs (MAC =0). For this reason it's a good idea to run min_mac(vcfR, min.mac=1)
#' before using a SNP dataset in downstream analyses.
#'
#' @param vcfR a vcfR object
#' @param min.mac an integer specifying the minimum minor allele count for a
#' SNP to be retained (e.g. 'min.mac=3' would remove all SNPs with a MAC of 2 or less)
#' @return if 'min.mac' is not specified, the allele frequency spectrum is returned.
#' If 'min.mac' is specified, SNPs falling below the MAC cutoff will be removed,
#' and the filtered vcfR object will be returned.
#' @examples
#' min_mac(vcfR=SNPfiltR::vcfR.example)
#' @export
min_mac <- function(vcfR, min.mac=NULL){

  #if vcfR is not class vcfR, fail gracefully
  if (class(vcfR) != "vcfR"){
    stop("specified vcfR object must be of class 'vcfR'")
  }

  if (is.null(min.mac)){

    #print message
    print("no filtering cutoff provided, vcf will be returned unfiltered")

    #convert vcfR to matrix and make numeric
    gt.matrix<-vcfR::extract.gt(vcfR)
    gt.matrix[gt.matrix == "0/0"]<-0
    gt.matrix[gt.matrix == "0/1"]<-1
    gt.matrix[gt.matrix == "1/1"]<-2
    class(gt.matrix) <- "numeric"

    #calc sfs
    sfs<-rowSums(gt.matrix, na.rm = TRUE)
    #fold sfs
    for (i in 1:length(sfs)) {
      if (sfs[i] <= sum(!is.na(gt.matrix[i,]))){}
      else {
        sfs[i]<-(sum(!is.na(gt.matrix[i,]))*2 - sfs[i])
      }
    }

    #hist folded mac with cutoff shown
    hist(sfs,
         main="folded SFS",
         xlab = "MAC")

    #return unfiltered vcf
    return(vcfR)

  }

  else{

    #convert vcfR to matrix and make numeric
    gt.matrix<-vcfR::extract.gt(vcfR)
    gt.matrix[gt.matrix == "0/0"]<-0
    gt.matrix[gt.matrix == "0/1"]<-1
    gt.matrix[gt.matrix == "1/1"]<-2
    class(gt.matrix) <- "numeric"

    #calc sfs
    sfs<-rowSums(gt.matrix, na.rm = TRUE)

    #fold sfs
    for (i in 1:length(sfs)) {
      if (sfs[i] <= sum(!is.na(gt.matrix[i,]))){}
      else {
        sfs[i]<-(sum(!is.na(gt.matrix[i,]))*2 - sfs[i])
      }
    }

    #hist folded mac with cutoff shown
    hist(sfs,
         main="folded SFS",
         xlab = "MAC")
    abline(v=min.mac-1,
           col="red")

    #calculate % of SNPs to be removed
    p<-round((sum(sfs < min.mac)/length(sfs))*100, 2)

    #print message to user
    print(paste0(p, "% of SNPs fell below a minor allele count of ", min.mac, " and were removed from the VCF"))

    #filter vcfR
    vcfR <- vcfR[sfs >= min.mac,]

    #return vcf object
    return(vcfR)

  }

}
