#bring in full vcf file
#read in vcf as vcfR
vcfR <- read.vcfR("~/Desktop/aph.data/unzipped.filtered.vcf")
dim(vcfR@gt)

#
head(vcfR@fix)
levels(as.factor(vcfR@fix[,1]))

distance.thin <- function(vcfR, distance=NULL){
#set distance
j=distance
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
  for (i in 2:length(fix.sub)){
    #store logical indicating whether this SNP is greater than j base pairs from the previous SNP
    k[i]<- fix[i] > prev+j
    #if it is, then we keep this SNP, making it the new 'previous' for assessing the next point.
    #If we dont keep the SNP, we don't update the closest point
    if (fix[i] > prev+j){
      prev<-fix[i]
    }
  }
  #now we subset the entire chromosome by the logical vector to retain only points > j bp apart
  keep<-fix.sub[k]
  #finally, we generate a logical vector indicating whether to keep each SNP while matching the original order of the chrom
  #we paste these into a continuously built vector for each chrom
  g<-c(g, fix.sub %in% keep)
}
#subset vcfR locus info based on logical
vcfR@fix<-vcfR@fix[g,]
#subset genotypes based on logical
vcfR@gt<-vcfR@gt[g,]
#return vcfR
return(vcfR)
}


