#' Hard filter a vcf file by depth and genotype quality (gq)
#'
#' This function requires a vcfR object as input.
#' There are slots for varying the n parameter across M-1, M, and M-1 (as recommended by Paris et al. 2017).
#' After running stacks with each of the n options, plug the output vcf files into this 
#' function to visualize the effect of varying m on number of SNPs and loci built to 
#' recognize which value optimizes the n parameter for your dataset at the 'R80' cutoff (Paris et al. 2017). 
#'
#' @param vcfR a vcfR object
#' @param depth an integer representing the minimum depth for genotype calls that you wish to retain
#' (e.g. 'depth = 5' would remove all genotypes with a sequencing depth of 4 reads or less)
#' @param gq an integer representing the minimum genotype quality for genotype calls that you wish to retain
#' #' (e.g. 'gq = 30' would remove all genotypes with a quality score of 29 or lower)
#' @return The vcfR object input, with the sites failing specified filters converted to 'NA'
#' @export
hard.filter.vcf <- function(vcfR, depth=NULL, gq=NULL){
  
  #extract depth from the vcf
  dp.matrix<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)
  
  #calculate the SNPs that fall below the depth filter
  i<-round((sum(dp.matrix < depth, na.rm = TRUE)/sum(!is.na(dp.matrix)))*100, 2)
  #report filter
  print(paste0(i,"% of genotypes fall below a read depth of ",depth," and were converted to NA"))
  
  #convert to NAs
  dp.matrix[dp.matrix < depth] <- NA
  vcfR@gt[,-1][ is.na(dp.matrix) == TRUE ] <- NA
  
  #extract gq from the vcf
  gq.matrix<- vcfR::extract.gt(vcfR, element='GQ', as.numeric=TRUE)
  
  #calculate the SNPs that fall below the gq filter
  j<-round((sum(gq.matrix < gq, na.rm = TRUE)/sum(!is.na(gq.matrix)))*100, 2)
  #report filter
  print(paste0(j,"% of genotypes fall below a genotype quality of ",gq," and were converted to NA"))
  
  #convert to NAs
  gq.matrix[gq.matrix < gq] <- NA
  vcfR@gt[,-1][ is.na(gq.matrix) == TRUE ] <- NA
  
  return(vcfR)
}

