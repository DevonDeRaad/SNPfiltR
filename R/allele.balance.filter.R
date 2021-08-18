#' Filter out heterozygous sites failing an allele balance check
#'
#' This function requires a vcfR object as input, and returns a vcfR object filtered 
#' to convert heterozygous sites with an allele balance falling outside of .25-.75 to 'NA'.
#' From: https://www.ddocent.com/filtering/ "Allele balance: a number between 0 and 1 
#' representing the ratio of reads showing the reference allele to all reads, 
#' considering only reads from individuals called as heterozygous, we expect that the 
#' allele balance in our data (for real loci) should be close to 0.5".
#'
#' @param vcfR a vcfR object
#' @return The vcfR object input, with the sites failing the allele balance filter converted to 'NA'
#' @export
filter.allele.balance <- function(vcfR){
  
  #extract AD from the vcf
  ad.matrix<- vcfR::extract.gt(vcfR, element='AD')
  #extract GT from the vcf
  gt.matrix<- vcfR::extract.gt(vcfR, element='GT')
  
  #mask dp matrix to include only called hets from gt matrix
  ad.matrix[gt.matrix != "0/1"]<-NA
  
  #split allele 1 depth from allele 2 depth
  al1<-structure(as.numeric(gsub(",.*", "", ad.matrix)), dim=dim(ad.matrix))
  al2<-structure(as.numeric(gsub(".*,", "", ad.matrix)), dim=dim(ad.matrix))
  
  #calculate percentage of hets failing AB filter
  AB<-al1/(al1 + al2) > .75 | al1/(al1 + al2) <.25
  p<-round(sum(AB, na.rm = TRUE) / sum(is.na(AB) == FALSE)*100, 2)
  j<-round(sum(AB, na.rm = TRUE) / sum(is.na(gt.matrix) == FALSE)*100, 2)
  
  print(paste0(p,"% of het genotypes (",j,"% of all genotypes) fall outside of .25 - .75 allele balance and were converted to NA"))
  
  #convert failing genotypes to NA
  vcfR@gt[,-1][AB]<-NA
  
  return(vcfR)
}
  
