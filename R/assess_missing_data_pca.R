#' Vizualise how missing data thresholds affect sample clustering
#'
#' This function can be run in two ways: 1) Without 'thresholds' specified. This will run a PCA
#' for the input vcf without filtering, and visualize the clustering of samples in two-dimensional
#' space, coloring each sample according to a priori population assignment given in the popmap.
#' 2) With 'thresholds' specified. This will filter your input vcf file to the specified
#' missing data thresholds, and run a PCA for each filtering iteration.
#' For each iteration, a 2D plot will be output showing clustering according to the
#' specified popmap. This option is ideal for assessing the effects of missing data
#' on clustering patterns.
#' @param vcfR a vcfR object
#' @param popmap set of population assignments that will be used to color code the plots
#' @param thresholds optionally specify a vector of missing data filtering thresholds to explore
#' @param clustering use partitioning around medoids (PAM) to do unsupervised
#' clustering on the output? (default = TRUE, max clusters = # of levels in popmap + 2)

#' @return a series of plots showing the clustering of all samples in two-dimensional space
#' @examples
#' assess_missing_data_pca(vcfR = SNPfiltR::vcfR.example,
#' popmap = SNPfiltR::popmap,
#' thresholds = c(.6,.8))
#' @export
assess_missing_data_pca <- function(vcfR,
                              popmap=NULL,
                              thresholds=NULL,
                              clustering=TRUE){

  #bind global variables
  PC1<- NULL
  PC2<- NULL
  pop<- NULL

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

      #record pam clustering info directly on PCA
      m=c()
      for (z in 2:(length(levels(as.factor(popmap$pop)))+2)){
        m[z]<-mean(cluster::silhouette(cluster::pam(pca.scores, z))[, "sil_width"])
      }

      #plot pam clustering info
      plot(m,
           type = "o",
           xlab = "K",
           ylab = "PAM silhouette",
           main=paste0("PCA PAM clustering"))

      #make dataframe
      pam.df<-data.frame(n.groups=2:(length(levels(as.factor(popmap$pop)))+2),
                         likelihood=m[-1])

      #run pam best clustering scheme
      pam.clust<-cluster::pam(pca.scores, pam.df$n.groups[pam.df$likelihood==max(pam.df$likelihood)])

      #match order for pop from popmap into this df
      pca.scores$pop<-popmap$pop[order(popmap$id == colnames(vcfR@gt)[-1])]
      pca.scores$pam.clust<-pam.clust$clustering
      #add missingness to df
      pca.scores$missing<-miss

      #record percentage of variance explained in the PCA
      var_frac <- pca$eig/sum(pca$eig)

      #plot pam clusters versus a priori clusters
      print(
        ggplot2::ggplot(pca.scores,
                        ggplot2::aes(x=PC1,
                                     y=PC2,
                                     color=pop,
                                     shape=as.factor(pam.clust))
        ) +
          ggplot2::ggtitle(paste0("PCA clustering analysis"))+
          ggplot2::xlab(paste0("PC1, ", round(var_frac[1]*100, 2), "% variance explained"))+
          ggplot2::ylab(paste0("PC2, ", round(var_frac[2]*100, 2), "% variance explained"))+
          ggplot2::geom_point(cex = 4,
                              alpha=.75)+
          ggplot2::theme_classic()+
          ggplot2::guides(color=ggplot2::guide_legend(title="popmap assignment"),
                          shape=ggplot2::guide_legend(title="PAM clusters"))
      )

      #plot PCA color coding by missing data percentage
      print(
        ggplot2::ggplot(pca.scores,
                        ggplot2::aes(x=PC1,
                                     y=PC2,
                                     color=missing)
        ) +
          ggplot2::ggtitle(paste0("PCA clustering analysis"))+
          ggplot2::geom_point(cex = 4,
                              alpha=.75)+
          ggplot2::labs(color = "proportion\nmissing data")+
          ggplot2::xlab(paste0("PC1, ", round(var_frac[1]*100, 2), "% variance explained"))+
          ggplot2::ylab(paste0("PC2, ", round(var_frac[2]*100, 2), "% variance explained"))+
          ggplot2::theme_classic()
      )

      #return df
      return(pca.scores)

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

        #record pam clustering info directly on PCA
        m=c()
        for (z in 2:(length(levels(as.factor(popmap$pop)))+2)){
          m[z]<-mean(cluster::silhouette(cluster::pam(pca.scores, z))[, "sil_width"])
        }

        #plot pam clustering info
        plot(m,
             type = "o",
             xlab = "K",
             ylab = "PAM silhouette",
             main=paste0(i*100,"% SNP completeness cutoff PAM clustering results"))

        #make dataframe
        pam.df<-data.frame(n.groups=2:(length(levels(as.factor(popmap$pop)))+2),
                           likelihood=m[-1])

        #run pam best clustering scheme
        pam.clust<-cluster::pam(pca.scores, pam.df$n.groups[pam.df$likelihood==max(pam.df$likelihood)])

        #match order for pop from popmap into this df
        pca.scores$pop<-popmap$pop[order(popmap$id == colnames(vcfR@gt)[-1])]
        pca.scores$pam.clust<-pam.clust$clustering
        #add missingness to df
        pca.scores$missing<-miss

        #record percentage of variance explained in the PCA
        var_frac <- pca$eig/sum(pca$eig)

        #plot pam clusters versus a priori clusters
        print(
          ggplot2::ggplot(pca.scores,
                          ggplot2::aes(x=PC1,
                                       y=PC2,
                                       color=pop,
                                       shape=as.factor(pam.clust))
          ) +
            ggplot2::ggtitle(paste0(i*100,"% SNP completeness cutoff PCA"))+
            ggplot2::geom_point(cex = 4,
                                alpha=.75)+
            ggplot2::theme_classic()+
            ggplot2::xlab(paste0("PC1, ", round(var_frac[1]*100, 2), "% variance explained"))+
            ggplot2::ylab(paste0("PC2, ", round(var_frac[2]*100, 2), "% variance explained"))+
            ggplot2::guides(color=ggplot2::guide_legend(title="popmap assignment"),
                            shape=ggplot2::guide_legend(title="PAM clusters"))
        )

        #plot PCA color coding by missing data percentage
        print(
          ggplot2::ggplot(pca.scores,
                          ggplot2::aes(x=PC1,
                                       y=PC2,
                                       color=missing)
          ) +
            ggplot2::ggtitle(paste0(i*100,"% SNP completeness cutoff PCA"))+
            ggplot2::geom_point(cex = 4,
                                alpha=.75)+
            ggplot2::labs(color = "proportion\nmissing data")+
            ggplot2::xlab(paste0("PC1, ", round(var_frac[1]*100, 2), "% variance explained"))+
            ggplot2::ylab(paste0("PC2, ", round(var_frac[2]*100, 2), "% variance explained"))+
            ggplot2::theme_classic()
        )

        #fill list with relevant dataframe
        dfs[[i]]<-pca.scores

        #clean objects
        pca.scores<-NULL
        var_frac<-NULL

        #close for loop
      }

    #return list
    return(dfs)

    #close else statement from if (thresholds==NULL)
  }

  } #close if statement from if (clustering == TRUE)

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

      #match order for pop from popmap into this df
      pca.scores$pop<-popmap$pop[order(popmap$id == colnames(vcfR@gt)[-1])]
      #add missingness to df
      pca.scores$missing<-miss

      #record percentage of variance explained in the PCA
      var_frac <- pca$eig/sum(pca$eig)

      #plot pam clusters versus a priori clusters
      print(
        ggplot2::ggplot(pca.scores,
                        ggplot2::aes(x=PC1,
                                     y=PC2,
                                     color=pop)
        ) +
          ggplot2::ggtitle(paste0("PCA clustering analysis"))+
          ggplot2::xlab(paste0("PC1, ", round(var_frac[1]*100, 2), "% variance explained"))+
          ggplot2::ylab(paste0("PC2, ", round(var_frac[2]*100, 2), "% variance explained"))+
          ggplot2::geom_point(cex = 4,
                              alpha=.75)+
          ggplot2::theme_classic()+
          ggplot2::guides(color=ggplot2::guide_legend(title="popmap assignment"))
      )

      #plot PCA color coding by missing data percentage
      print(
        ggplot2::ggplot(pca.scores,
                        ggplot2::aes(x=PC1,
                                     y=PC2,
                                     color=missing)
        ) +
          ggplot2::ggtitle("PCA clustering anlysis")+
          ggplot2::geom_point(cex = 4,
                              alpha=.75)+
          ggplot2::labs(color = "proportion\nmissing data")+
          ggplot2::xlab(paste0("PC1, ", round(var_frac[1]*100, 2), "% variance explained"))+
          ggplot2::ylab(paste0("PC2, ", round(var_frac[2]*100, 2), "% variance explained"))+
          ggplot2::theme_classic()
      )

      #return df
      return(pca.scores)

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

        #match order for pop from popmap into this df
        pca.scores$pop<-popmap$pop[order(popmap$id == colnames(vcfR@gt)[-1])]
        #add missingness to df
        pca.scores$missing<-miss

        #record percentage of variance explained in the PCA
        var_frac <- pca$eig/sum(pca$eig)

        #plot pam clusters versus a priori clusters
        print(
          ggplot2::ggplot(pca.scores,
                          ggplot2::aes(x=PC1,
                                       y=PC2,
                                       color=pop)
          ) +
            ggplot2::ggtitle(paste0(i*100,"% SNP completeness cutoff PCA"))+
            ggplot2::geom_point(cex = 4,
                                alpha=.75)+
            ggplot2::theme_classic()+
            ggplot2::xlab(paste0("PC1, ", round(var_frac[1]*100, 2), "% variance explained"))+
            ggplot2::ylab(paste0("PC2, ", round(var_frac[2]*100, 2), "% variance explained"))+
            ggplot2::guides(color=ggplot2::guide_legend(title="popmap assignment"))
        )

        #plot PCA color coding by missing data percentage
        print(
          ggplot2::ggplot(pca.scores,
                          ggplot2::aes(x=PC1,
                                       y=PC2,
                                       color=missing)
          ) +
            ggplot2::ggtitle(paste0(i*100,"% SNP completeness cutoff PCA"))+
            ggplot2::geom_point(cex = 4,
                                alpha=.75)+
            ggplot2::labs(color = "proportion\nmissing data")+
            ggplot2::xlab(paste0("PC1, ", round(var_frac[1]*100, 2), "% variance explained"))+
            ggplot2::ylab(paste0("PC2, ", round(var_frac[2]*100, 2), "% variance explained"))+
            ggplot2::theme_classic()
        )

        #fill list with relevant dataframe
        dfs[[i]]<-pca.scores

        #clean objects
        pca.scores<-NULL
        var_frac<-NULL

        #close for loop
      }

      #return list
      return(dfs)

      #close else statement from if (thresholds==NULL)
    }

  } #close if statement from if (clustering == FALSE)

  #close function
}
