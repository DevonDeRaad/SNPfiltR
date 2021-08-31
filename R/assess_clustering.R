#' Vizualise how missing data thresholds affect sample clustering
#'
#' This function can be run in two ways: 1) Without 'min.mac' specified. This will run DAPC for a
#' minimum MAC of 1,2,3,4,5,10, visualize the results, and tell you how the DAPC clustering
#' compares to your a priori population assignments, and how many SNPs are retained at each MAC cutoff.
#' 2) With 'min.mac' specified. This will print your folded Site Frequency Spectrum (SFS) and
#' show you where your specified min. MAC count falls. It will then return your vcfR object with SNPs
#' falling below your min. MAC threshold removed.
#' MAF cutoffs can be helpful in removing spurious and uninformative loci from the dataset, but also
#' have the potential to bias downstream inferences. Linck and Battey (2019) have an excellent paper on
#' just this topic. From the paper- "We recommend researchers using model‐based programs to describe
#' population structure observe the following best practices: (a) duplicate analyses with nonparametric
#' methods suchas PCA and DAPC with cross validation (b) exclude singletons (c) compare alignments with
#' multiple assembly parameters When seeking to exclude only singletons in alignments with missing data
#' (a ubiquitous problem for reduced‐representation library preparation methods), it is preferable to filter
#' by the count (rather than frequency) of the minor allele, because variation in the amount of missing data
#' across an alignment will cause a static frequency cutoff to remove different SFS classes at different sites".
#' Based on the differences between DAPC runs and the variation in number of SNPs retained, you should
#' be able to make an informed decision about whether to implement a MAC cutoff for your dataset.
#' Note: previous filtering steps may have converted called genotypes to 'NAs' and resulting in invariant SNPs
#' (MAC =0) for this reason it's a good idea to at least run min.mac(vcfR, min.mac=1).
#' @param vcfR a vcfR object
#' @param thresholds a vector specifying the filtering thresholds you want to explore
#' @param popmap set of population assignments that will be used to color code the plots
#' @return a series of plots showing the clustering of all samples in two-dimensional space
#' @examples
#' assess_clustering(
#' vcfR=system.file("extdata","unfiltered.vcf.gz",package="SNPfiltR",mustWork=TRUE),
#' thresholds=c(60,80,100))
#' @export
assess_clustering <- function(vcfR, popmap=NULL, thresholds=NULL){

  #if vcfR is not class vcfR, fail gracefully
  if (class(vcfR) != "vcfR"){
    stop("specified vcfR object must be of class 'vcfR'")
  }

  if (is.null(popmap)){
    stop("popmap must be provided in order to compare clustering of a priori defined groups")
  }

  if (class(popmap) != "data.frame"){
    stop("popmap must be of class 'data.frame'")
  }

    if (colnames(popmap)[1] != "id" | colnames(popmap)[2] != "pop"){
      stop("popmap must be a dataframe with two columns, 'id' and 'pop'")
  }

  if (is.null(thresholds)) {

    #run clustering with no filters

    #convert vcfR into genlight
    genlight<-vcfR::vcfR2genlight(vcfR)

    #run dapc for mac 1,2,3,4,5,10
    for (i in c(1,2,3,4,5,10)){

      #filter genlight by given mac
      genlight<-genlight[,sfs >= i]
      #subset sfs vector to only samples left in the vcf
      sfs<-sfs[sfs >= i]

      #assign samples to the number of groups present in popmap, retain all PCAs
      grp<-adegenet::find.clusters(genlight, n.pca = ncol(gt.matrix)-1, n.clust = length(levels(popmap$pop)))

      #check how well that assignment matched up to the provided popmap
      samps<-merge(popmap, data.frame(group=grp$grp, id=labels(grp$grp)), by='id')
      print(paste0("for ", i, " minimum MAC cutoff, compare k means clustering to popmap assignment"))
      print(table(samps$pop, samps$group))

      #run dapc, retain all discriminant axes, and enough PC axes to explain 75% of variance
      dapc1<-adegenet::dapc(genlight, grp$grp, n.da = length(levels(popmap$pop))-1, pca.select = "percVar", perc.pca = 75)

      #plot compoplot
      adegenet::compoplot(dapc1, legend=FALSE, show.lab =TRUE, cex.names=.4, main=paste0("min. MAC ",i,", total SNPs ",length(sfs)))

      #print
      print(paste0("DAPC with min. MAC ", i, " and ", length(sfs), " total SNPs, complete"))

      return(vcfR)

    }

  }
  else {

    #run clustering at each specified filtering level

    return(vcfR)

  }

}
