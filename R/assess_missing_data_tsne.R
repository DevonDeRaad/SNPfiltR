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
#' @param popmap set of population assignments that will be used to color code the plots
#' @param thresholds a vector specifying the missing data filtering thresholds to explore
#' @param perplexity numerical value specifying the perplexity paramter during t-SNE (default: 10)
#' @param iterations a numerical value specifying the number of iterations for t-SNE
#' (default: 1000)
#' @param initial_dims a numerical value specifying the number of
#' initial_dimensions for t-SNE (default: 5)
#'  @param clustering use partitioning around medoids (PAM) to do unsupervised
#' clustering on the output? (default = TRUE, max clusters = # of levels in popmap + 2)

#' @return a series of plots showing the clustering of all samples in two-dimensional space
#' @examples
#' assess_clustering_tsne(
#' vcfR=system.file("extdata","unfiltered.vcf.gz",package="SNPfiltR",mustWork=TRUE),
#' thresholds=c(60,80,100))
#' @export
assess_missing_data_tsne <- function(vcfR,
                              popmap=NULL,
                              thresholds=NULL,
                              perplexity=NULL,
                              iterations=NULL,
                              initial_dims=NULL,
                              clustering=TRUE){

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

  #ensure popmap naming conventions followed
  if (colnames(popmap)[1] != "id"){
    stop("popmap must be a dataframe with column 1 named 'id'")
  }

  #ensure popmap naming conventions followed
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

  #set perplexity to default (10), or user specified value
  if (is.null(perplexity)){
    w<-10
  } else{
    w<-perplexity
  }

  #set iterations to default (1000), or user specified value
  if (is.null(iterations)){
    q<-1000
  } else{
    q<-iterations
  }

  #set initial_dims to default (5), or user specified value
  if (is.null(initial_dims)){
    p<-5
  } else{
    p<-initial_dims
  }

  #this PCA cannot tolerate invariant SNP positions, so check for invariant SNP positions
  #convert vcfR to matrix and make numeric
  gt.matrix<-vcfR::extract.gt(vcfR)
  gt.matrix[gt.matrix == "0/0"]<-0
  gt.matrix[gt.matrix == "0/1"]<-1
  gt.matrix[gt.matrix == "1/1"]<-2
  class(gt.matrix) <- "numeric"

  #calc sfs
  sfs<-rowSums(gt.matrix, na.rm = TRUE)

  #calc number of invariant SNPs
  g<-sum(sfs < 1)

  #If there are invariant SNPs, kill the function, and tell user that invariant SNPs aren't allowed
  if (g != 0){
    stop("invariant SNPs detected in input vcf. Invariant sites must be filtered prior to input")
  }

  #calculate missingness by individual
  miss<-colSums(is.na(gt.matrix))/nrow(gt.matrix)

  ###
  ###
  ###
  #clustering module
  ###
  ###
  ###

  #if clustering = TRUE, start here:
  if (clustering == TRUE){

  #if checks on inputs pass, and clustering == TRUE, and thresholds are not specified, start here:
    if (is.null(thresholds)) {

    #run clustering with no filters
    #convert vcfR into genlight
    gen<-vcfR::vcfR2genlight(vcfR)

    #execute PCA using this genlight
    #retain number of PC axes equivalent to the number of populations being discriminated + 2
    pca<-adegenet::glPca(gen,
                         nf=length(levels(as.factor(popmap$pop)))+2)

    #pull pca scores out of df
    pca.scores<-as.data.frame(pca$scores)
    #pca.scores$pop<-popmap$pop

    ###############################################
    ###############################################
    # t-SNE
    ###############################################
    ###############################################

    # t-SNE on PCA
    tsne_p5<- tsne::tsne(pca.scores,
                         max_iter=q,
                         perplexity=w,
                         initial_dims=p)

    #save as dataframe
    tsne.df<-as.data.frame(tsne_p5)

    #record pam clustering info
    m=c()
    for (z in 2:(length(levels(as.factor(popmap$pop)))+2)){
      m[z]<-mean(cluster::silhouette(cluster::pam(tsne_p5, z))[, "sil_width"])
    }

    #plot pam clustering info
    plot(m,
         type = "o",
         xlab = "K",
         ylab = "PAM silhouette",
         main=paste0("t-SNE PAM clustering"))

    #make dataframe
    pam.df<-data.frame(n.groups=2:(length(levels(as.factor(popmap$pop)))+2),
                       likelihood=m[-1])

    #run pam best clustering scheme
    pam.clust<-cluster::pam(tsne_p5, pam.df$n.groups[pam.df$likelihood==max(pam.df$likelihood)])

    #match order for pop from popmap into this df
    tsne.df$pop<-popmap$pop[order(popmap$id == colnames(vcfR@gt)[-1])]
    tsne.df$pam.clust<-pam.clust$clustering
    #add missingness to df
    tsne.df$missing<-miss

    #plot pam clusters versus a priori clusters
    print(
      ggplot2::ggplot(tsne.df,
                      ggplot2::aes(x=V1,
                                   y=V2,
                                   color=pop,
                                   shape=as.factor(pam.clust))
      ) +
        ggplot2::ggtitle(paste0("t-SNE clustering analysis"))+
        ggplot2::geom_point(cex = 4,
                            alpha=.75)+
        ggplot2::theme_classic()
    )

    #plot t-SNE color coding by missing data percentage
    print(
      ggplot2::ggplot(tsne.df,
                      ggplot2::aes(x=V1,
                                   y=V2,
                                   color=missing)
      ) +
        ggplot2::ggtitle(paste0("t-SNE clustering analysis"))+
        ggplot2::geom_point(cex = 4,
                            alpha=.75)+
        ggplot2::labs(color = "proportion\nmissing data")+
        ggplot2::theme_classic()
    )

    #make clean df with info
    df<-data.frame(tsne.ax1=tsne.df$V1,
                   tsne.ax2=tsne.df$V2,
                   popmap.pop=popmap$pop,
                   pam.pop=as.factor(pam.clust$clustering),
                   missing=tsne.df$missing)

    #return
    return(df)

    #close if(is.null(thresholds))
  }

  #if thresholds are provided, run loop over all filtering thresholds here:
  else{

    #open list to fill meta info
    dfs<-list()

    for (i in thresholds){

      #if specified cutoff is not between 0-1 fail gracefully
      if (i < 0 | i > 1){
        stop("specified threshold must be a proportion between 0 and 1")
      }

      #otherwise filter vcfR based on given threshold
      vcfR.filt<-SNPfiltR::missing_by_snp(vcfR = vcfR,
                                          cutoff = i)

      #convert vcfR into genlight
      gen<-vcfR::vcfR2genlight(vcfR.filt)

      #execute PCA using this genlight
      #retain number of PC axes equivalent to the number of populations being discriminated + 2
      pca<-adegenet::glPca(gen,
                           nf=length(levels(as.factor(popmap$pop)))+2)

      #pull pca scores out of df
      pca.scores<-as.data.frame(pca$scores)

      ###############################################
      ###############################################
      # t-SNE
      ###############################################
      ###############################################

      # t-SNE on PCA output
      tsne_p5<- tsne::tsne(pca.scores,
                           max_iter=q,
                           perplexity=w,
                           initial_dims = p)

      #save as data frame
      tsne.df<-as.data.frame(tsne_p5)

      #record pam clustering info
      m=c()
      for (z in 2:(length(levels(as.factor(popmap$pop)))+2)){
        m[z]<-mean(cluster::silhouette(cluster::pam(tsne_p5, z))[, "sil_width"])
      }

      #plot pam clustering info
      plot(m,
           type = "o",
           xlab = "K",
           ylab = "PAM silhouette",
           main=paste0("t-SNE ",i*100,"% SNP completeness cutoff PAM clustering"))

      #make dataframe
      pam.df<-data.frame(n.groups=2:(length(levels(as.factor(popmap$pop)))+2),
                         likelihood=m[-1])

      #run pam best clustering scheme
      pam.clust<-cluster::pam(tsne_p5, pam.df$n.groups[pam.df$likelihood==max(pam.df$likelihood)])

      #match order for pop from popmap into this df
      tsne.df$pop<-popmap$pop[order(popmap$id == colnames(vcfR.filt@gt)[-1])]
      tsne.df$pam.clust<-pam.clust$clustering
      #add missingness to df
      tsne.df$missing<-miss

      #plot pam clusters versus a priori clusters
      print(
        ggplot2::ggplot(tsne.df,
                        ggplot2::aes(x=V1,
                                     y=V2,
                                     color=pop,
                                     shape=as.factor(pam.clust))
        ) +
          ggplot2::ggtitle(paste0("t-SNE clustering analysis ",i*100,"% SNP completeness cutoff"))+
          ggplot2::geom_point(cex = 4,
                              alpha=.75)+
          ggplot2::theme_classic()
      )

      #plot t-SNE color coding by missing data percentage
      print(
        ggplot2::ggplot(tsne.df,
                        ggplot2::aes(x=V1,
                                     y=V2,
                                     color=missing)
        ) +
          ggplot2::ggtitle(paste0("t-SNE clustering analysis ",i*100,"% SNP completeness cutoff"))+
          ggplot2::geom_point(cex = 4,
                              alpha=.75)+
          ggplot2::labs(color = "proportion\nmissing data")+
          ggplot2::theme_classic()
      )

      #fill list with relevant dataframe
      dfs[[i]]<-tsne.df

      #clean df object
      tsne.df<-NULL

    #close for loop
    }
    return(dfs)
  #close else statement
  }

  #close if(clustering==TRUE){} statement
  }

  ###
  ###
  ###
  #non-clustering module
  ###
  ###
  ###

  #If clustering has been turned off, just make PCAs colored by popmap
  if (clustering == FALSE){

    #if checks on inputs pass, and clustering == TRUE, and thresholds are not specified, start here:
    if (is.null(thresholds)) {

      #convert vcfR into genlight
      gen<-vcfR::vcfR2genlight(vcfR)

      #execute PCA using this genlight
      #retain number of PC axes equivalent to the number of populations being discriminated + 2
      pca<-adegenet::glPca(gen,
                           nf=length(levels(as.factor(popmap$pop)))+2)

      #pull pca scores out of df
      pca.scores<-as.data.frame(pca$scores)
      #pca.scores$pop<-popmap$pop

      ###############################################
      ###############################################
      # t-SNE
      ###############################################
      ###############################################

      # t-SNE on PCA
      tsne_p5<- tsne::tsne(pca.scores,
                           max_iter=q,
                           perplexity=w,
                           initial_dims=p)

      #save as dataframe
      tsne.df<-as.data.frame(tsne_p5)

      #match order for pop from popmap into this df
      tsne.df$pop<-popmap$pop[order(popmap$id == colnames(vcfR@gt)[-1])]

      #add missingness to df
      tsne.df$missing<-miss

      #plot pam clusters versus a priori clusters
      print(
        ggplot2::ggplot(tsne.df,
                        ggplot2::aes(x=V1,
                                     y=V2,
                                     color=pop)
        ) +
          ggplot2::ggtitle(paste0("t-SNE clustering analysis"))+
          ggplot2::geom_point(cex = 4,
                              alpha=.75)+
          ggplot2::theme_classic()
      )

      #plot t-SNE color coding by missing data percentage
      print(
        ggplot2::ggplot(tsne.df,
                        ggplot2::aes(x=V1,
                                     y=V2,
                                     color=missing)
        ) +
          ggplot2::ggtitle(paste0("t-SNE clustering analysis"))+
          ggplot2::geom_point(cex = 4,
                              alpha=.75)+
          ggplot2::labs(color = "proportion\nmissing data")+
          ggplot2::theme_classic()
      )

      #make clean df with info
      df<-data.frame(tsne.ax1=tsne.df$V1,
                     tsne.ax2=tsne.df$V2,
                     popmap.pop=popmap$pop)

      #return
      return(df)

      #close if(is.null(thresholds))
    }

    #if thresholds are provided, run loop over all filtering thresholds here:
    else{

      #open list to fill meta info
      dfs<-list()

      for (i in thresholds){

        #if specified cutoff is not between 0-1 fail gracefully
        if (i < 0 | i > 1){
          stop("specified threshold must be a proportion between 0 and 1")
        }

        #otherwise filter vcfR based on given threshold
        vcfR.filt<-SNPfiltR::missing_by_snp(vcfR = vcfR,
                                            cutoff = i)

        #convert vcfR into genlight
        gen<-vcfR::vcfR2genlight(vcfR.filt)

        #execute PCA using this genlight
        #retain number of PC axes equivalent to the number of populations being discriminated + 2
        pca<-adegenet::glPca(gen,
                             nf=length(levels(as.factor(popmap$pop)))+2)

        #pull pca scores out of df
        pca.scores<-as.data.frame(pca$scores)

        ###############################################
        ###############################################
        # t-SNE
        ###############################################
        ###############################################

        # t-SNE on PCA output
        tsne_p5<- tsne::tsne(pca.scores,
                             max_iter=q,
                             perplexity=w,
                             initial_dims = p)

        #save as data frame
        tsne.df<-as.data.frame(tsne_p5)

        #match order for pop from popmap into this df
        tsne.df$pop<-popmap$pop[order(popmap$id == colnames(vcfR.filt@gt)[-1])]

        #add missingness to df
        tsne.df$missing<-miss

        #plot pam clusters versus a priori clusters
        print(
          ggplot2::ggplot(tsne.df,
                          ggplot2::aes(x=V1,
                                       y=V2,
                                       color=pop)
          ) +
            ggplot2::ggtitle(paste0("t-SNE clustering analysis ",i*100,"% SNP completeness cutoff"))+
            ggplot2::geom_point(cex = 4,
                                alpha=.75)+
            ggplot2::theme_classic()
        )

        #plot t-SNE color coding by missing data percentage
        print(
          ggplot2::ggplot(tsne.df,
                          ggplot2::aes(x=V1,
                                       y=V2,
                                       color=missing)
          ) +
            ggplot2::ggtitle(paste0("t-SNE clustering analysis ",i*100,"% SNP completeness cutoff"))+
            ggplot2::geom_point(cex = 4,
                                alpha=.75)+
            ggplot2::labs(color = "proportion\nmissing data")+
            ggplot2::theme_classic()
        )

        #fill list with relevant dataframe
        dfs[[i]]<-tsne.df

        #clean df object
        tsne.df<-NULL

        #close for loop
      }
      return(dfs)
      #close else statement
    }

    #close if(clustering==FALSE){} statement
  }

  #close function
}
