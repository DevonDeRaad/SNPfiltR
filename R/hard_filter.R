#' Hard filter a vcf file by depth and genotype quality (gq)
#'
#' This function requires a vcfR object as input.
#' The user can then specify the minimum value for depth of coverage required to retain
#' a called genotype (must be numeric).
#' Additionally, the user can specify a minimum genotype quality required to retain
#' a called genotype (again, must be numeric).
#'
#' @param vcfR a vcfR object
#' @param depth an integer representing the minimum depth for genotype calls that you
#' wish to retain
#' (e.g. 'depth = 5' would remove all genotypes with a sequencing depth of 4 reads or less)
#' @param gq an integer representing the minimum genotype quality for
#' genotype calls that you wish to retain
#' (e.g. 'gq = 30' would remove all genotypes with a quality score of 29 or lower)
#' @return The vcfR object input, with the sites failing specified filters converted to 'NA'
#' @examples
#' hard_filter(vcfR = SNPfiltR:::vcfR.example, depth = 5)
#' hard_filter(vcfR = SNPfiltR:::vcfR.example, depth = 5, gq = 30)
#' @export
hard_filter <- function(vcfR, depth=NULL, gq=NULL){

  #if specified vcfR is not class 'vcfR', fail gracefully
  if (class(vcfR) != "vcfR"){
    stop("specified vcfR object must be of class 'vcfR'")
  }

  #if depth is specified, start here
  if (!is.null(depth)) {

    if (is.numeric(depth) != "TRUE"){
      stop("specified depth cutoff must be numeric")
    }

  #extract depth from the vcf
  dp.matrix<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)

  #calculate the SNPs that fall below the depth filter
  i<-round((sum(dp.matrix < depth, na.rm = TRUE)/sum(!is.na(dp.matrix)))*100, 2)
  #report filter
  print(paste0(i,"% of genotypes fall below a read depth of ",depth," and were converted to NA"))

  #convert to NAs
  dp.matrix[dp.matrix < depth] <- NA
  vcfR@gt[,-1][ is.na(dp.matrix) == TRUE ] <- NA
  #close if statement
  }

  #if no depth is specified
  else{
    #print user message
    print("no depth cutoff provided, exploratory visualization will be generated.")

    #extract depth from the vcf
    dp.matrix<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)

    #set plotting parameters
    #plot histogram of depth
    hist(dp.matrix,
         xlab = "genotype depth")
    abline(v=mean(dp.matrix, na.rm = TRUE),
           col="red",
           lty="dashed")

    #zoomed in histogram
    hist(dp.matrix[dp.matrix < 25],
         xlab = "genotype depth")

  }

  #if gq is specified
  if (!is.null(gq)) {

    if (is.numeric(gq) != "TRUE"){
      stop("specified genotype quality cutoff must be numeric")
    }

  #extract gq from the vcf
  gq.matrix<- vcfR::extract.gt(vcfR, element='GQ', as.numeric=TRUE)

  #calculate the SNPs that fall below the gq filter
  j<-round((sum(gq.matrix < gq, na.rm = TRUE)/sum(!is.na(gq.matrix)))*100, 2)

  #report filter
  print(paste0(j,"% of genotypes fall below a genotype quality of ",gq," and were converted to NA"))

  #convert to NAs
  gq.matrix[gq.matrix < gq] <- NA
  vcfR@gt[,-1][ is.na(gq.matrix) == TRUE ] <- NA

  #close if statement
  }

  else{
    print("no genotype quality cutoff provided, exploratory visualization will be generated.")

    #extract gq from the vcf
    gq.matrix<- vcfR::extract.gt(vcfR, element='GQ', as.numeric=TRUE)

    #set plotting parameters
    par(mfrow=c(2,1))
    #plot histogram of depth
    hist(gq.matrix,
         xlab = "genotype quality")
    abline(v=mean(gq.matrix, na.rm = TRUE),
           col="red",
           lty="dashed")

  }

  #return
  return(vcfR)

#close function
}

