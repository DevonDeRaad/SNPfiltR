#' Filter a vcf file based on distance between SNPs on a given scaffold
#'
#' This function requires a vcfR object as input, and returns a vcfR object filtered
#' to retain only SNPs greater than a specified distance apart on each scaffold.
#' The function starts by automatically retaining the first SNP on a given scaffold,
#' and then subsequently keeping the next SNP that is greater than the specified distance away,
#' until it reaches the end of the scaffold/chromosome. This function scales well with
#' an increasing number of SNPs, but poorly with an increasing number of scaffolds/chromosomes.
#' For this reason, there is a built in progress bar,
#' to monitor potentially long-running executuions with many scaffolds.
#' This type of filtering is often employed to reduce linkage among input SNPs.
#'
#' @param vcfR a vcfR object
#' @return An identical vcfR object, except that SNPs separated by less than the
#' specified distance have been removed from the file
#' @examples
#' distance_thin(vcfR = SNPfiltR:::vcfR.example, min.distance = 1000)
#' @export
distance_thin <- function(vcfR, min.distance=NULL){

  #if vcfR is not class vcfR, fail gracefully
  if (class(vcfR) != "vcfR"){
    stop("specified vcfR object must be of class 'vcfR'")
  }

  #function is useless if no distance is provided
  if (is.null(min.distance)){
    stop("filtering distance must be provided")
  }

  #logical test specifying that the minimum distance between SNPs for filtering must be at least 1 base pair or the logic doesn't work
  if (min.distance < 1){
    stop("filtering distance must be >= 1 base pair")
  }

  #logical test to ensure formatting of input vcf will allow accurate analysis of position in genome
  if (colnames(vcfR@fix)[1] != "CHROM"){
    stop("vcfR incorrectly formatted. vcfR@fix column 1 must be 'CHROM'")
  }

  #logical test to ensure formatting of input vcf will allow accurate analysis of position in genome
  if (colnames(vcfR@fix)[2] != "POS"){
    stop("vcfR incorrectly formatted. vcfR@fix column 2 must be 'POS'")
  }

#set min distance specified by user
j=min.distance

#generate dataframe containing information for chromosome and bp locality of each SNP
df<-as.data.frame(vcfR@fix[,1:2])

#generate list of all unique chromosomes in alphabetical order
chroms<-levels(as.factor(df$CHROM))

#intialize empty df to hold filtering
keep.df<-data.frame()

#make progress bar
pb <- txtProgressBar(min = 0, max = length(chroms), style = 3)

#begin tracker
pbtrack<-1

#loop over each chromosome
#for t in vector containing the name of each chromosome
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

      #else, use a for loop to determine which SNPs to keep that satisfy our distance criteria
      for (i in 2:length(fix.sub)){

        #store logical indicating whether this SNP is greater than j base pairs from the previous SNP
        k[i]<- fix.sub[i] > prev+j

        #if it is, then we keep this SNP, making it the new 'previous' for assessing the next point.
        #If we don't keep the SNP, we don't update the closest point
          if (fix.sub[i] > prev+j){
          prev<-fix.sub[i]
          }

      #close for loop
      }
    #close else statement
    }

  #make a dataframe with the precise info for each SNP for this chromosome
  chrom.df<-data.frame(CHROM=rep(t, times=length(fix.sub)), POS=fix.sub, keep=k)

  #now we rbind in the information for this chromosome to the overall df
  keep.df<-rbind(keep.df,chrom.df)

  #empty df for this chrom to prepare for the next one
  chrom.df<-NULL

  #update progress bar
  setTxtProgressBar(pb, pbtrack)

  #update tracker
  pbtrack<-pbtrack+1

} #close for loop, start over on next chromosome

#close progress bar
close(pb)

#order the dataframe to match the order of the input vcf file
order.df<-keep.df[match(paste(df$CHROM,df$POS), paste(keep.df$CHROM,keep.df$POS)),]

#subset vcfR locus info based on the logical column from our dataframe
#vcfR@fix<-vcfR@fix[order.df$keep,]
#subset genotypes based on logical
#vcfR@gt<-vcfR@gt[order.df$keep,]

#realized there is no need to do this subsetting separately
vcfR<-vcfR[order.df$keep,]

#calculate number of total SNPs input
z<-nrow(keep.df)

#calculate total SNPs retained
p<-nrow(vcfR@fix)

#print info to screen
print(paste0(p," out of ",z," input SNPs were not located within ",j," base-pairs of another SNP and were retained despite filtering"))

#return vcfR
return(vcfR)
}


