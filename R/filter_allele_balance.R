#' Filter out heterozygous genotypes failing an allele balance check
#'
#' This function requires a vcfR object as input, and returns a vcfR object filtered
#' to convert heterozygous sites with an allele balance falling outside of the
#' specified ratio to 'NA'. If no ratio is specified, a default .25-.75 limit is imposed.
#' From the [dDocent filtering page](https://www.ddocent.com/filtering/) "Allele balance:
#' a number between 0 and 1 representing the ratio of reads showing the reference allele to
#' all reads, considering only reads from individuals called as heterozygous, we expect that the
#' allele balance in our data (for real loci) should be close to 0.5".
#'
#' @param vcfR a vcfR object
#' @param min.ratio minumum allele ratio for a called het
#' @param max.ratio maximum allele ratio for a called het
#' @return An identical vcfR object, except that genotypes failing the allele balance
#' filter have been converted to 'NA'.
#' @examples
#' filter_allele_balance(vcfR = SNPfiltR::vcfR.example)
#' @export
filter_allele_balance <- function(vcfR,
                                  min.ratio=NULL,
                                  max.ratio=NULL){

  #if specified vcfR is not class 'vcfR', fail gracefully
  if (class(vcfR) != "vcfR"){
    stop("specified vcfR object must be of class 'vcfR'")
  }

  #set default parameter for minimum allele ratio for a het call (.25)
  if(is.null(min.ratio)){
    min.ratio=.25
  }

  #set default parameter for maximum allele ratio for a het call (.75)
  if(is.null(max.ratio)){
    max.ratio=.75
  }
  #if these ratios have been specified by the user, leave them alone

  #extract allele depth from the vcf

  #if allele depth is specified as 'AD', extract matrix
  ad.matrix<- vcfR::extract.gt(vcfR, element='AD')

  #write a test to identify how allele depth is specified, or not specified
  if (length(grep("AD",vcfR@gt[,1])) > 0.5){
    #if allele depth is specified as 'AD', extract matrix
    ad.matrix<- vcfR::extract.gt(vcfR, element='AD')
  }else if(length(grep("CATG",vcfR@gt[,1])) > 0.5){
    #print warning that this will be slow
    print("Warning, allele depth encoded only as raw CATG counts, must first index out relevant allele depths, which is time consuming")
    #if allele depth is specified as 'CATG', extract matrix of 'CATG' values
    full.matrix<- vcfR::extract.gt(vcfR, element='CATG')
    #open empty matrix to hold relevant AD info same size as full.matrix
    ad.matrix<-matrix(, nrow = nrow(full.matrix), ncol = ncol(full.matrix))
    #extract only REF and ALT allele depth values
    for (i in 1:nrow(full.matrix)){
        #convert reference allele to indexing tool, C > 1, A > 2, T > 3, G > 4
        k<-as.numeric(gsub("G","4",gsub("T","3",gsub("A","2",gsub("C","1",as.data.frame(vcfR@fix)$REF[i])))))
        #convert ALT allele to indexing tool, C > 1, A > 2, T > 3, G > 4
        l<-as.numeric(gsub("G","4",gsub("T","3",gsub("A","2",gsub("C","1",as.data.frame(vcfR@fix)$ALT[i])))))
        #extract the REF and ALT depths for each genotype and save them as a single vector, w/ values separated by a comma. Write that vector to the new matrix
        ad.matrix[i,]<-sapply(lapply(strsplit(full.matrix[i,],","), '[', c(k,l)), paste0, collapse=",")
    }
    print("CATG format converted to REF,ALT allele depth")
  }else{
    stop("allele depth is not specified in input vcf in 'AD' or 'CATG' format, therefore allele balance cannot be calculated")
  }

  #extract GT from the vcf
  gt.matrix<- vcfR::extract.gt(vcfR, element='GT')

  #mask ad matrix to include only called hets from gt matrix
  ad.matrix[gt.matrix != "0/1" & gt.matrix != "1/0"]<-NA

  #split allele 1 depth from allele 2 depth
  al1<-structure(as.numeric(gsub(",.*", "", ad.matrix)), dim=dim(ad.matrix))
  al2<-structure(as.numeric(gsub(".*,", "", ad.matrix)), dim=dim(ad.matrix))

  #calculate AB for each sample
  al.bal<-al1/(al1 + al2)

  #calculate logical storing whether each het genotype passes the filter
  AB<-al1/(al1 + al2) > max.ratio | al1/(al1 + al2) < min.ratio

  #calculate percent of het genotypes failing the filter
  p<-round(sum(AB, na.rm = TRUE) / sum(is.na(AB) == FALSE)*100, 2)

  #calculate overall percentage of genotypes failing the filter
  j<-round(sum(AB, na.rm = TRUE) / sum(is.na(gt.matrix) == FALSE)*100, 2)

  #print to user
  message(p,"% of het genotypes (",j,"% of all genotypes) fall outside of ",min.ratio," - ",max.ratio, " allele balance ratio and were converted to NA")

  #convert failing genotypes to NA
  vcfR@gt[,-1][AB]<-NA

  #make histogram of allele balance at all het genotypes
  graphics::hist(al.bal,
                 xlim = c(0,1),
                 ylab = "number of genotypes",
                 xlab = "Allele balance",
                 main ="allele balance distribution")
  graphics::abline(v=c(min.ratio,max.ratio),
                   col="red")

  #return vcfR
  return(vcfR)

  #close function
}
