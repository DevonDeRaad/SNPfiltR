#' Vizualise missing data per sample, remove samples above a missing data cutoff
#'
#' This function can be run in two ways: 1) Without 'cutoff' specified. This will vizualise the
#' amount of missing data in each sample across a variety of potential missing data cutoffs.
#' Additionally, it will show you a dotplot ordering the amount of overall missing data in each sample.
#' Based on these visualizations, you can make an informed decision on what you think might be an optimal
#' cutoff to remove samples that are missing too much data to be retained for downstream analyses. 2)
#' with 'cutoff' specified. This option will show you the dotplot with the cutoff you set, and then
#' remove samples above the missing data cutoff you set, and return the filtered vcf to you.
#'
#' Note: This decision is highly project specific, but these visualizations should help you get a feel for how
#' very low data samples cannot be rescued simply by a missing data SNP filter. If you want to remove
#' specific samples from your vcf that cannot be specified with a simple cutoff refer to this great tutorial:
#' https://knausb.github.io/vcfR_documentation/sequence_coverage.html which was the basis for the code
#' implemented in this function.
#'
#' @param vcfR a vcfR object
#' @param popmap if specifies, it must be a two column dataframe with columns names 'id' and 'pop'.
#' IDs must match the IDs in the vcfR object
#' @param cutoff a numeric value between 0-1 specifying the maximum proportion of missing data
#' allowed in a sample to be retained for downstream analyses
#' @return if 'cutoff' is not specified, will return a dataframe containing the average depth and proportion
#' missing data in each sample. If 'cutoff' is specified, the samples falling above the missing data cutoff
#' will be removed, and the filtered vcfR object will be returned.
#' @export
missing.by.sample <- function(vcfR, popmap=NULL, cutoff=NULL){

  #extract depth from the vcf
  dp<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)

        if (is.null(cutoff)){

            if (is.null(popmap)) {
                print("No popmap provided")
                }
            else {
              if (colnames(popmap)[1] != "id" | colnames(popmap)[2] != "pop"){
                stop("popmap must be a dataframe with two columns, 'id' and 'pop'")
              }
                #calculate missingness by pop here and make dotplot
                  #popmap must be a two column dataframe with 'id' and 'pop' columns
                  #id's must match the ids in the vcf file
                  #calculate missingness by individual
                  miss<-colSums(is.na(dp))/nrow(dp)
                  #calculate avg depth by individual
                  avg.depth<-colMeans(dp, na.rm = TRUE)
                  #store ordered column names
                  samples<-colnames(dp)
                  #create df
                  df.x<-data.frame(id=samples,missingness=miss,avg.depth=avg.depth, row.names = NULL)
                  #bring in popmap
                  df.f<-merge(df.x, popmap, by="id")
                  #plot missingness and depth by pop
                  plot1<-ggplot2::ggplot(df.f, ggplot2::aes(x= reorder(pop, -missingness), y=missingness)) +
                    ggplot2::geom_violin()+
                    ggplot2::theme_classic()+
                    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))+
                    ggplot2::ylab("proportion missing")+
                    ggplot2::geom_dotplot(binaxis='y', stackdir='center', dotsize=.8)
                  #dotplot avg depth of sequencing
                  plot2<-ggplot2::ggplot(df.f, ggplot2::aes(x= reorder(pop, -missingness), y=avg.depth)) +
                    ggplot2::geom_violin()+
                    ggplot2::theme_classic()+
                    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))+
                    ggplot2::ylab("average depth")+
                    ggplot2::geom_dotplot(binaxis='y', stackdir='center', dotsize=.8)

                  #
                  gridExtra::grid.arrange(plot1,plot2)

                }


    #initialize dataframe
    df.y<- data.frame(indiv=character(), snps.retained=numeric(), filt=numeric())
    #loop
    for (i in c(.5,.6,.7,.8,.9,1)){
      #get vector of individuals we are looping over
      indiv<-colnames(dp)
      #calculate the completeness cutoff for each snp to be retained in this iteration
      filt<-rep(i, times =ncol(dp))
      #calculate the number of loci successfully genotyped in each sample in this iteration
      snps.retained<-colSums(is.na(dp[(rowSums(is.na(dp))/ncol(dp) <= 1-i),]) == "FALSE")
      #append the data from this iteration to existing df
      df.y<-rbind(df.y, as.data.frame(cbind(indiv, filt, snps.retained)))
      #close for loop
    }

    #make columns numeric for plotting
    df.y$filt<-as.numeric(as.character(df.y$filt))
    df.y$snps.retained<-as.numeric(as.character(df.y$snps.retained))
    rownames(df.y)<-NULL
    #visualize color coded by individual
    print(
      ggplot2::ggplot(df.y, ggplot2::aes(x=filt, y=snps.retained, color=indiv))+
        ggplot2::geom_point()+
        ggplot2::ggtitle("SNPs retained by filtering scheme") +
        ggplot2::xlab("fraction of non-missing genotypes required to retain each SNP (0-1)")+
        ggplot2::ylab("non-missing SNPs retained per sample")+
        ggplot2::theme_light()+
        ggplot2::theme(legend.position="none")
    )

    #calculate missingness by individual
    miss<-colSums(is.na(dp))/nrow(dp)
    #show plot with missingness by sample
    dotchart(sort(miss), cex=.5, xlab = "proportion missing data")
    abline(v=c(.5,.6,.7,.8,.9,1), lty="dashed")

    #return the dataframe showing the missingness and avg depth per individual
    return(df.y)

      }
  else {

    #calculate missingness by individual
    miss<-colSums(is.na(dp))/nrow(dp)
    #vis plot to show where cutoff was set
    dotchart(sort(miss), cex=.5)
    abline(v=cutoff, col="red")

    #print
    print(paste0(length(labels(miss)[miss > cutoff])," samples are above a ",cutoff," missing data cutoff, and were removed from VCF"))
    #drop those individuals from vcfR
    vcfR@gt <- vcfR@gt[,!(colnames(vcfR@gt) %in% labels(miss)[miss > cutoff])]

    return(vcfR)
  }

}

