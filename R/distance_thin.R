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
#' distance_thin(vcfR = SNPfiltR:::vcfR.example, min.distance = 1000)
#' @export
distance_thin <- function(vcfR, min.distance=NULL){
  if (is.null(min.distance)){
    stop("filtering distance must be provided")
  }

  if (min.distance < 1){
    stop("filtering distance must be >= 1")
  }

  if (colnames(vcfR@fix[,1:2]) != c("CHROM", "POS")){
    stop("vcfR incorrectly formatted. vcfR@fix columns 1 & 2 must be 'CHROM' and 'POS'")
  }

#set distance
j=min.distance
#generate df of chrom and bp
df<-as.data.frame(vcfR@fix[,1:2])
#generate list of all unique chromosomes in alphabetical order
chroms<-levels(as.factor(df$CHROM))
#intialize empty df to hold filtering
keep.df<-data.frame()

#loop over each chromosome
#for t in 1:length of the number of unique chromosomes
for (t in chroms){
  #isolate the SNP positions on the given chromosome
  fix.sub<-as.numeric(df$POS[df$CHROM == t])
  #order the positions numerically
  fix.sub<-fix.sub[order(fix.sub)]
  #set the first position
  prev<-fix.sub[1]
  #initialize empty vector
  k<-c()
  #always keep first SNP on the chromosome
  k[1]<-TRUE
  #loop to decide whether to keep each following SNP
    if (length(fix.sub) < 2){
      #if chrom only has 1 SNP, do nothing
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

  #make a dataframe with the precise info for each SNP for this chromosome
  chrom.df<-data.frame(CHROM=rep(t, times=length(fix.sub)), POS=fix.sub, keep=k)

  #now we rbind in the information for this chromosome to the overall df
  keep.df<-rbind(keep.df,chrom.df)

  #empty df for this chrom to prepare for the next one
  chrom.df<-NULL

} #close for loop, start over on next chromosome

#order the dataframe to match the order of the input vcf file
order.df<-keep.df[match(paste(df$CHROM,df$POS), paste(keep.df$CHROM,keep.df$POS)),]

#subset vcfR locus info based on the logical column from our dataframe
vcfR@fix<-vcfR@fix[order.df$keep,]
#subset genotypes based on logical
vcfR@gt<-vcfR@gt[order.df$keep,]

#calculate number of total SNPs input
z<-nrow(keep.df)

#calculate total SNPs retained
p<-nrow(vcfR@fix)

print(paste0(p," out of ",z," input SNPs were not located within ",j," base-pairs of another SNP and were retained despite filtering"))

#return vcfR
return(vcfR)
}


