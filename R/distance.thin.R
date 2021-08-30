#' Filter a vcf file based on distance between SNPs on a given scaffold
#'
#' This function requires a vcfR object as input, and returns a vcfR object filtered
#' to retain only SNPs greater than a specified distance apart on each scaffold.
#' This type of filtering is often employed to reduce linkage amongst SNPs.
#'
#' @param vcfR a vcfR object
#' @return An identical vcfR object, except that SNPs separated by less than the
#' specified distance have been removed from the file
#' @examples
#' distance.thin(vcfR = SNPfiltR:::vcfR.example, min.distance = 1000)
#' @export
distance.thin <- function(vcfR, min.distance=NULL){
  if (is.null(min.distance)){
    stop("filtering distance must be provided")
  }

  if (min.distance < 1){
    stop("filtering distance must be >= 1")
  }

#set distance
j=min.distance
#generate all bp positions as a numeric vector
fix<-as.numeric(vcfR@fix[,2])
#intialize empty vector to hold filtering
g<-c()

#loop over each chromosome
for (t in 1:length(levels(as.factor(vcfR@fix[,1])))){
  #isolate the given chromosome
  fix.sub<-fix[vcfR@fix[,1] == levels(as.factor(vcfR@fix[,1]))[t]]
  #order the positions
  fix.sub<-fix.sub[order(fix.sub)]
  #set the first position
  prev<-fix.sub[1]
  #initialize empty vector
  k<-c()
  #always keep first SNP on the chromosome
  k[1]<-TRUE
  #loop to decide whether to keep each following SNP
    if (length(fix.sub) < 2){ #if chrom only has 1 SNP, do nothing
    } else{
      for (i in 2:length(fix.sub)){
        #store logical indicating whether this SNP is greater than j base pairs from the previous SNP
        k[i]<- fix.sub[i] > prev+j
        #if it is, then we keep this SNP, making it the new 'previous' for assessing the next point.
        #If we don't keep the SNP, we don't update the closest point
          if (fix.sub[i] > prev+j){
          prev<-fix.sub[i]
          } #close if statement
      } #close for loop
    } #close else statement

  #now we subset the entire chromosome by the logical vector to retain only points > j bp apart
  keep<-fix.sub[k]
  #finally, we generate a logical vector indicating whether to keep each SNP while matching the original order of the chrom
  #we paste these into a continuously built vector for each chrom
  g<-c(g, fix.sub %in% keep)
} #close for loop, start over on next chromosome

#calculate total SNPs input
z<-length(g)

#subset vcfR locus info based on logical
vcfR@fix<-vcfR@fix[g,]
#subset genotypes based on logical
vcfR@gt<-vcfR@gt[g,]

#calculate total SNPs retained
p<-nrow(vcfR@fix)

print(paste0(p," out of ",z," input SNPs were not located within ",j," base-pairs of another SNP and were retained despite filtering"))

#return vcfR
return(vcfR)
}


