

#' Generate Bimodality Index (BI) for microarray type of data.   
#'
#' The bimodality index is a continuous measure of the extent to 
#'  which a set of (univariate) data fits a two-component mixture model.
#'  The score is larger if the two components are balanced in size or if
#'  the separation between the two modes is larger.
#'  
#' @param dataset A matrix or data.frame, usually with columns representing samples
#'  and rows representing genes or proteins.
#' @param verbose A logical value; should the function output a stream of information
#'  while it is working?
#'  
#' @return a data frame containing six columns, with the rows corresponding to the
#'  rows of the original data set.  The columns contain the four parameters from the
#'  normal mixture model (mu1, mu2, sigma, and pi) along with the standardized
#'  distance 'delta' and the bimodal index (BI).
#'  
#' @seealso \code{\link{bimodalIndex}} which this function wraps.
#' @export
#' @importFrom mclust mclustBIC
#' @examples
#' data(lung.dataset)
#' bi <- bimodalIndex(lung.dataset, verbose=FALSE)
#' summary(bi)
bimodalIndex <- function(dataset, verbose=TRUE) {
  bim <- matrix(NA, nrow=nrow(dataset), ncol=6)
  if (verbose) cat("1 ")
  for (i in 1:nrow(dataset)) {
    if (verbose && 0 == i%%100) cat(".")
    if (verbose && 0 == i%%1000) cat(paste("\n", 1 + i/1000, ' ', sep=''))
    x <- as.vector(as.matrix(dataset[i, ]))
    if (any(is.na(x))) next
    mc <- mclust::Mclust(x, G = 2, modelNames = "E")
    sigma <- sqrt(mc$parameters$variance$sigmasq)
    delta <- abs(diff(mc$parameters$mean))/sigma
    #    pi <- max(mc$parameters$pro)
    pi <- mc$parameters$pro[1]
    bi <- delta * sqrt(pi*(1-pi))
    bim[i,] <-  c(mc$parameters$mean, sigma=sigma, delta=delta, pi=pi, bim=bi)
  }
  if(verbose) cat("\n")
  dimnames(bim) <- list(rownames(dataset),
                        c("mu1", "mu2", "sigma", "delta", "pi", "BI"))
  bim <- as.data.frame(bim)
  return(bim)
}


#' Generate density-plot with Bimodality Index (BI) for microarray type of data.   
#'
#' In processing microarry data, sometime we want to find some genes with meaningful 
#'  expression values that can separate into at least two groups. Here we wrap the
#'  function 'bimodalIndex' in ClassDiscovery and function 'Mclust' in mclust and plot the
#'  density with the grouping info.
#'  
#' @param x input a numeric vector that used to calculate bimodality index (BI).
#' @param G input how many groups users want to separate based on mclust.
#' @param annCol input an information of samples that is used to combine the sample. 
#' @param annColors input the specified colors of the sample information.
#' @param xlab character string specifying the label of X axis.
#' @param main character string specifying the graph title.
#'  
#' @return a density-plot of BI value in the title of the legend with the grouping info.
#'  
#' @seealso \code{\link{bimodalIndex}} which this function wraps.
#' @export
#' @importFrom mclust mclustBIC
#' @examples
#' data(lung.dataset)
#' BI_densityplot(x=lung.dataset["AFFX-r2-Ec-bioD-5_at",])
#' BI_densityplot(x=lung.dataset["AFFX-r2-Ec-bioD-5_at",],
#'                annCol=c("1"="AFFX-Low", "2"="AFFX-High"),
#'                xlab="AFFX-r2-Ec-bioD-5_at")
#' BI_densityplot(x=lung.dataset["AFFX-r2-Ec-bioD-5_at",],
#'                annCol=c("1"="AFFX-Low", "2"="AFFX-High"),
#'                annColors=c("AFFX-Low"="blue", "AFFX-High"="red"))
#' BI_densityplot(x=lung.dataset[2,], G=3, annCol=c("1"="Low", "2"="Mid", "3"="High"))                
BI_densityplot <- function(x, G=2, annCol=c("1"="Low", "2"="High"), annColors=NULL,
                           xlab="Expression Value (log2)",
                           main="BI-Density_Plot"){
  #require(mclust)
  if (G!=length(annCol)) {
    stop("Number of Groups (G) must have the same length of annotation (annCol).\n")
  }
  if (G==length(annCol)) {
    bi <- bimodalIndex(as.matrix(t(x)), verbose=FALSE)$BI
    bscore <- plyr::revalue(as.character(mclust::Mclust(x, G=G, "E")$classification),
                            annCol)
    d_long <- data.frame(x, bscore)
    density.x <- min(x)
    density.y <- max(density(x)$y)
  
    if (!is.null(annColors)){
      p1 <- ggplot2::ggplot(d_long, ggplot2::aes(x=x)) +
        ggplot2::geom_density() + 
        ggplot2::geom_rug(data=d_long, ggplot2::aes(x=x, y=0, colour=bscore)) +
        ggplot2::scale_color_manual(values = annColors)
    }
    if (is.null(annColors)){
      p1 <- ggplot2::ggplot(d_long, ggplot2::aes(x=x)) +
        ggplot2::geom_density() +
        ggplot2::geom_rug(data=d_long, ggplot2::aes(x=x, y=0, colour=bscore))
    }
    p1 <- p1 + ggplot2::theme_bw() +
      ggplot2::theme(legend.position="right") +
      ggplot2::xlab(xlab) +  ggplot2::ylab("Density") +
      ggplot2::ggtitle(main) +
      ggplot2::labs(colour=paste("BI:", round(bi, 3)))
    print(p1)
  }
}


#' Separate the groups based on Bimodality Index (BI) for microarray type of data.   
#'
#' In processing microarry data, sometime we want to find some genes with meaningful 
#'  expression values that can separate into at least two groups. Here we wrap the
#'  function 'bimodalIndex' in ClassDiscovery and function 'Mclust' in mclust and return
#'  the the groups.
#'  
#' @param x input a numeric vector that used to calculate bimodality index (BI).
#' @param G input how many groups users want to separate based on mclust.
#' @param annCol input an information of samples that is used to combine the sample.
#'   
#' @return a character vector for the groups from 'mclust'.
#'  
#' @seealso \code{\link{Mclust}} which this function wraps.
#' @export
#' @importFrom mclust mclustBIC
#' @examples
#' data(lung.dataset)
#' BI_groups(x=lung.dataset["AFFX-r2-Ec-bioD-5_at",])
#' BI_groups(x=lung.dataset["AFFX-r2-Ec-bioD-5_at",],
#'            G=2, annCol=c("1"="AFFX-Low", "2"="AFFX-High"))
#' BI_groups(x=lung.dataset[2,], G=3, annCol=c("1"="Low", "2"="Mid", "3"="High"))           
BI_groups <- function(x, G=2, annCol=c("1"="Low", "2"="High")){
  #require(mclust)
  if (G!=length(annCol)) {
    stop("Number of Groups (G) must have the same length of annotation (annCol).\n")
  }
  if (G==length(annCol)) {
    bscore <- plyr::revalue(as.character(mclust::Mclust(x, G=G, "E")$classification),
                            annCol)
  }
  return(bscore)
}






