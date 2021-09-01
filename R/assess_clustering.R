#' Vizualise how missing data thresholds affect sample clustering
#'
#' This function can be run in two ways: 1) Without 'thresholds' specified. This will run t-SNE
#' for the input vcf without filtering, and visualize the clustering of samples in two-dimensional
#' space, coloring each sample according to a priori population assignment given in the popmap.
#' 2) With 'thresholds' specified. This will filter your input vcf file to the specified
#' missing data thresholds, and run a t-SNE clustering analysis for each filtering iteration.
#' For each iteration, a 2D plot will be output showing clustering according to the
#' specified popmap. This option is ideal for assessing the effects of missing data
#' on clustering patterns.
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

  #if popmap not provided, fail gracefully
  if (is.null(popmap)){
    stop("popmap must be provided in order to compare clustering of a priori defined groups")
  }

  #if popmap not formatted correctly, fail gracefully
  if (class(popmap) != "data.frame"){
    stop("popmap must be of class 'data.frame'")
  }

  #ensure naming conventions followed
  if (colnames(popmap)[1] != "id"){
    stop("popmap must be a dataframe with column 1 named 'id'")
  }

  #ensure naming conventions followed
  if (colnames(popmap)[2] != "pop"){
    stop("popmap must be a dataframe with column 2 named 'pop'")
  }

  #check that id column length in popmap matches the number of samples in the vcf file
  if (length(popmap$id) != length(colnames(vcfR@gt))-1){
    stop("popmap ID's must match exactly the ID's in input vcf")
  }

  #check that id's match the ids in the vcf file
  if (all(popmap$id %in% colnames(vcfR@gt)) == FALSE){
    stop("popmap ID's must match exactly the ID's in input vcf")
  }

  #if checks pass, and thresholds not provided, start here:
  if (is.null(thresholds)) {

    #run clustering with no filters

    #convert vcfR into genlight
    genlight<-vcfR::vcfR2genlight(vcfR)

    ###############################################
    ###############################################
    # t-SNE
    ###############################################
    ###############################################

    # prepare plot labels and such
    # this makes it so it is grouped by DAPC clusters
    colors = rainbow(length(unique(results$grp)))
    names(colors) = unique(results$grp)
    ecb = function(x,y){plot(x,t='n'); text(x, labels=results$grp, col=colors[results$grp])}

    # t-SNE on principal components of scaled data
    # adjust perplexity, initial_dims
    # can do k=3 for 3D plot
    # should do only <50 variables
    # can do it on pca$li (if you reduce the number of components), or on cmdsplot2$points
    tsne_p5 = tsne(pca1$tab, epoch_callback=ecb, max_iter=5000, perplexity=5, initial_dims=5)

    # tSNE plot with DAPC groups
    plot(tsne_p5, main="t-SNE perplexity=5 with DAPC optimal k and clusters", col=results$grp, pch=16)

    # pam clustering with optimal k from DAPC
    for (i in 2:10){
      print(paste(i, mean(silhouette(pam(tsne_p5, i))[, "sil_width"])))
    }
    #pam prefers same 6 groups as DAPC
    pam(tsne_p5, 6)

    # determine optimal k of tSNE via hierarchical clustering with BIC
    # adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
    tsne_p5_clust <- Mclust(tsne_p5)
    mclust_grps_tsne_p5 <- as.numeric(tsne_p5_clust$classification)
    max(mclust_grps_tsne_p5)
    # t-SNE p5 with optimal k and clusters of RF via hierarchical clustering
    plot(tsne_p5,
         xlab="Scaling Dimension 1",
         ylab="Scaling Dimension 2",
         main="t-SNE p5 RF optimal K and clusters (hierarchical clustering)",
         col=mclust_grps_tsne_p5, pch=16)
    mclust_grps_tsne_p5
    f<-as.data.frame(tsne_p5)
    # tSNE with optimal k and clusters via hierarchical clustering
    ggplot(data=f, aes(x=V1,
                       y=V2,
                       col=as.factor(mclust_grps_tsne_p5)))+
      geom_point(cex=3)+
      theme_classic()

    cbind(rownames(pca1$tab), mclust_grps_tsne_p5)
      #return vcfR
      return(vcfR)

      #close if statement
    }

    #if thresholds are provided, run loop over all filtering thresholds here:
  }
  else {

    #run clustering at each specified filtering level

    return(vcfR)

  }

}
