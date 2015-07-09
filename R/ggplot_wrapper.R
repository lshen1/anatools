#' Fit a beta-uniform mixture model to a set of p-values.  
#'
#' In ClassComparison package in R there is a `Bum` class which is used to fit
#'  a beta-uniform mixture model to a set of p-values. This function wraps several
#'  functions, including 'Bum', 'hist', 'cutoffSignificant', and 'selectSignificant'
#'  into one function and presents the results in one figure.
#' @param pvals numeric vector containing values between '0' and '1'.
#' @param method either the false discovery rate (if by = 'FDR') or the posterior
#'  probability (if by = 'EmpiricalBayes'). Default is set to 'FDR'. 
#' @param fdrs input the false discovery rate (FDR).
#' @param res positive integer scalar specifying the resolution at which to plot
#'  the fitted distribution curve.
#' @param heights specify the heights of Bum plot and FDR table.
#' @param main character string specifying the graph title.    
#' @return a figure with BUM plot on the top and the FDR table on the bottom.
#' @seealso \code{\link{Bum}}, \code{\link{hist}}, \code{\link{cutoffSignificant}},
#'  and \code{\link{selectSignificant}} which this function wraps.
#' @export
#' @examples
#' fake.data <- c(runif(700), rbeta(300, 0.3, 1))
#' gg_Bum_Tab(fake.data, main="", heights=c(0.75, 0.25))
gg_Bum_Tab <- function(pvals, method='FDR', fdrs=c(0.01, 0.05, 0.1, 0.2),
                       res=100, heights=c(0.7, 0.3), main="Bum Plot") {
  a <- ClassComparison::Bum(pvals)
  counts <- sapply(fdrs, function(alpha) ClassComparison::countSignificant(a, alpha, by=method))
  cuts <- sapply(fdrs, function(alpha) ClassComparison::cutoffSignificant(a, alpha, by=method))
  mat <- data.frame(fdrs, counts, cuts)
  colnames(mat) <- c("FDR", "Number Significant", "P-value Cutoff")
  
  fdat <- data.frame(pvals=pvals)
  xvals <- (0:res)/res
  fit <- a@lhat + (1 - a@lhat) * dbeta(xvals, a@ahat, 1)
  betaf <- data.frame(xvals=xvals, fit=fit)[-1,]
  
  den.plot <- ggplot2::ggplot(fdat, aes(x=pvals)) + 
    geom_histogram(aes(y = ..density..), colour="black", fill = "white", binwidth = 0.01) +
    geom_line(data=betaf, aes(xvals,fit), colour="darkgreen") +
    geom_hline(yintercept = a@pihat, colour="blue") +
    labs(x="P Values", y="Density", title=main) +
    theme_classic()
  
  tab.plot <- gridExtra::tableGrob(mat)
  
  gridExtra::grid.arrange(den.plot, tab.plot, ncol=1, heights=heights, main="")
}


#' Generate a plot of table shows the features whose p-values less than cutoffs.  
#'
#' In ClassComparison package in R there is a `Bum` class which is used to fit
#'  a beta-uniform mixture model to a set of p-values. However sometimes the number
#'  of features are limitted, we would like to expand the feature list by ignoring
#'  FDR correstion.
#' @param pvals numeric vector containing values between '0' and '1'.
#' @param annot input the name of features used in the figure.
#' @param cutoffs input the cutoff of the p-values.
#' @param main character string specifying the graph title.    
#' @return a figure of table shows the features whose p-values less than cutoffs.
#' @seealso \code{\link{Bum}}, \code{\link{hist}}, \code{\link{cutoffSignificant}},
#'  and \code{\link{selectSignificant}} which this function wraps.
#' @export
#' @examples
#' fake.data <- c(rbeta(100, 0.8, 1))
#' annot.name <- paste("Probe", 1:100, sep="_")
#' gg_pval_Tab(fake.data, annot.name)
#' gg_pval_Tab(fake.data, annot.name, main="Tumor vs. Normal", footnote="Without FDR correction.")
gg_pval_Tab <- function(pvals, annot, cutoffs=c(0.01, 0.05, 0.1),
                        main="", footnote="") {
  counts <- sapply(cutoffs, function(alpha) length(which(pvals <= alpha))) 
  items <- sapply(cutoffs, function(alpha) paste(annot[which(pvals <= alpha)], collapse=", ")) 
  
  mat <- data.frame(cutoffs, counts, items)
  colnames(mat) <- c("p-values", "Number Significant", "Features")
  
  ### rearrange layout of the contents in Feature:
  d = sapply(lapply(mat$Features, strwrap, width=40), paste, collapse="\n")
  mat$Features <- d
  
  ### start new plot:
  grid.newpage()
  grid::pushViewport(viewport(height=0.8,width=0.95))
  
  ### plot table:
  p <- gridExtra::tableGrob(mat)
  
  ### plot title, footnote, and roundrect:
  h <- grobHeight(p)
  w <- grobWidth(p)
  title <- textGrob(main, y=unit(0.5,"npc") + 0.4*h, 
                    vjust=0, gp=gpar(fontsize=20))
  footnote <- textGrob(footnote, 
                       x=unit(0.5,"npc") - 0.5*w,
                       y=unit(0.5,"npc") - 0.5*h, 
                       vjust=1, hjust=0,gp=gpar( fontface="italic"))
  rect <- roundrectGrob(y=unit(0,"line"), height=unit(1,"npc") +unit(1,"lines"),
                        just="bottom", r=unit(0.05, "snpc"))
  gt <- gTree(children=gList(rect, p, title, footnote))
  
  ### generate plot:
  grid.draw(gt)
}



#' Generate a plot of table shows the features whose coefficient greater than cutoffs.  
#'
#' We like to list the features whose coefficient greater than certain cutoffs in a table.
#' @param cefs numeric vector containing values between '0' and '1'.
#' @param annot input the name of features used in the figure.
#' @param cutoffs input the cutoff of the p-values.
#' @param main character string specifying the graph title.    
#' @return a figure of table shows the features whose p-values less than cutoffs.
#' @seealso \code{\link{Bum}}, \code{\link{hist}}, \code{\link{cutoffSignificant}},
#'  and \code{\link{selectSignificant}} which this function wraps.
#' @export
#' @examples
#' fake.data <- c(rbeta(10, 0.8, 1))
#' annot.name <- paste("Probe", 1:100, sep="_")
#' gg_cef_Tab(fake.data, annot.name)
#' gg_cef_Tab(fake.data, annot.name, main="Tumor vs. Normal", footnote="Without FDR correction.")
gg_cef_Tab <- function(cefs, annot, cutoffs=c(0.5, 0.6, 0.7, 0.8),
                        main="", footnote="") {
  counts <- sapply(cutoffs, function(alpha) length(which(cefs >= alpha))) 
  items <- sapply(cutoffs, function(alpha) paste(annot[which(cefs >= alpha)], collapse=", ")) 
  
  mat <- data.frame(cutoffs, counts, items)
  colnames(mat) <- c("Coefficient", "Number Significant", "Features")
  
  ### rearrange layout of the contents in Feature:
  d = sapply(lapply(mat$Features, strwrap, width=40), paste, collapse="\n")
  mat$Features <- d
  
  ### start new plot:
  grid.newpage()
  grid::pushViewport(viewport(height=0.8,width=0.95))
  
  ### plot table:
  p <- gridExtra::tableGrob(mat)
  
  ### plot title, footnote, and roundrect:
  h <- grobHeight(p)
  w <- grobWidth(p)
  title <- textGrob(main, y=unit(0.5,"npc") + 0.4*h, 
                    vjust=0, gp=gpar(fontsize=20))
  footnote <- textGrob(footnote, 
                       x=unit(0.5,"npc") - 0.5*w,
                       y=unit(0.5,"npc") - 0.5*h, 
                       vjust=1, hjust=0,gp=gpar( fontface="italic"))
  rect <- roundrectGrob(y=unit(0,"line"), height=unit(1,"npc") +unit(1,"lines"),
                        just="bottom", r=unit(0.05, "snpc"))
  gt <- gTree(children=gList(rect, p, title, footnote))
  
  ### generate plot:
  grid.draw(gt)
}






#' Generate a plot of table shows the results of Fisher's Exact Test for Count Data.  
#'
#' Performs Fisher's exact test for testing the null of independence of rows and columns
#'  in a contingency table with fixed marginals.
#' @param x either a two-dimensional contingency table in matrix form, or a factor object.
#' @param y a factor object; ignored if 'x' is a matrix.
#' @param workspace an integer specifying the size of the workspace used in the network algorithm.
#'   In units of 4 bytes.  Only used for non-simulated p-values larger than 2 by 2 tables.
#' @param hybrid a logical. Only used for larger than 2 by 2 tables, in which cases it indicates
#'  whether the exact probabilities (default) or a hybrid approximation thereof should be computed.
#'   See 'Details'.
#' @param control a list with named components for low level algorithm control. At present the only
#'  one used is '"mult"', a positive integer >= 2 with default 30 used only for larger than 2 by 2
#'  tables. This says how many times as much space should be allocated to paths as to keys: see
#'  file 'fexact.c' in the sources of this package.
#' @param or the hypothesized odds ratio.  Only used in the 2 by 2 case.
#' @param alternative indicates the alternative hypothesis and must be one of '"two.sided"',
#'  '"greater"' or '"less"'.  You can specify just the initial letter.  Only used in the 2 by 2 case.
#' @param conf.int logical indicating if a confidence interval for the odds ratio in a 2 by 2 table
#'  should be computed (and returned).
#' @param conf.level confidence level for the returned confidence interval. Only used in the 2 by 2
#'  case and if 'conf.int = TRUE'.
#' @param simulate.p.value a logical indicating whether to compute p-values by Monte Carlo simulation,
#'  in larger than 2 by 2 tables.
#' @param B: an integer specifying the number of replicates used in the Monte Carlo test. 
#' @return a figure of table shows the results of Fisher's Exact Test for Count Data.
#' @seealso \code{\link{fisher.test}} which this function wraps.
#' @export
#' @examples
#'TeaTasting <- matrix(c(3, 1, 1, 3),
#'                     nrow = 2, dimnames = list(Guess = c("Milk", "Tea"),
#'                     Truth = c("Milk", "Tea")))
#'TeaTasting                     
#'fisher.test(TeaTasting)
#'gg_fisher.test_Tab(TeaTasting, show.dimnames=TRUE)
#'Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
#'              dimnames = list(income = c("< 15k", "15-25k", "25-40k", "> 40k"),
#'                              satisfaction = c("VeryD", "LittleD", "ModerateS", "VeryS"))
#'              )
#'Job
#'fisher.test(Job)
#'gg_fisher.test_Tab(Job)
gg_fisher.test_Tab <- function(x, y = NULL,
                               main="Fisher's Exact Test for Count Data",
                               show.dimnames=FALSE,
                               ...) {
  ### parsing x or x/y:
  DNAME <- deparse(substitute(x))
  METHOD <- "Fisher's Exact Test for Count Data"
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (is.matrix(x)) {
    if (any(dim(x) < 2L)) 
      stop("'x' must have at least 2 rows and columns")
    if (!is.numeric(x) || any(x < 0) || anyNA(x)) 
      stop("all entries of 'x' must be nonnegative and finite")
    if (!is.integer(x)) {
      xo <- x
      x <- round(x)
      if (any(x > .Machine$integer.max)) 
        stop("'x' has entries too large to be integer")
      if (!identical(TRUE, (ax <- all.equal(xo, x)))) 
        warning(gettextf("'x' has been rounded to integer: %s", 
                         ax), domain = NA)
      storage.mode(x) <- "integer"
    }
  }
  else {
    if (is.null(y)) 
      stop("if 'x' is not a matrix, 'y' must be given")
    if (length(x) != length(y)) 
      stop("'x' and 'y' must have the same length")
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) 
      stop("'x' and 'y' must have at least 2 levels")
    x <- table(x, y)
  }
  
  if(show.dimnames==TRUE){
    rownames(x) <- paste(names(attributes(x)$dimnames)[1], attributes(x)$dimnames[[1]], sep=".")
    colnames(x) <- paste(names(attributes(x)$dimnames)[2], attributes(x)$dimnames[[2]], sep=".")
  }
  
  ### perform fisher.test:
  res <- stats::fisher.test(x, ...)
  footnote <- paste("Data: ", DNAME, "\n",
                    "p-value = ", round(res$p.value, digits=4), "\n",
                    "Alternative hypothesis: ", res$alternative, "\n",
                    ifelse(is.null(res$conf.int), paste("", "\n", sep=""),
                           paste("95 percent confidence interval:", "\n",
                                 round(res$conf.int[1], digits=4), " - ",
                                 round(res$conf.int[2], digits=4), "\n", sep="")),
                    ifelse(is.null(res$estimate), paste("", "\n", sep=""),
                           paste("odds ratio:", "\n", round(res$estimate, digits=4), sep="")),
                    sep="")
  
  ### start new plot:
  grid.newpage()
  grid::pushViewport(viewport(height=0.8,width=0.95))
  
  ### plot table:
  p <- gridExtra::tableGrob(x)
  
  ### plot title, footnote, and roundrect:
  h <- grobHeight(p)
  w <- grobWidth(p)
  title <- textGrob(main, y=unit(0.5,"npc") + 0.4*h, 
                    vjust=0, gp=gpar(fontsize=20))
  footnote <- textGrob(footnote, 
                       x=unit(0.5,"npc") - 0.5*w,
                       y=unit(0.5,"npc") - 0.5*h, 
                       vjust=1, hjust=0,gp=gpar( fontface="italic"))
  rect <- roundrectGrob(y=unit(0,"line"), height=unit(1,"npc") +unit(1,"lines"),
                        just="bottom", r=unit(0.05, "snpc"))
  gt <- gTree(children=gList(rect, p, title, footnote))
  
  ### generate plot:
  grid.draw(gt)
 
}


#' Generate density-plot for QC accessment on microarray type of data. 
#'
#' In processing microarry data, usually we need to see how the expression values on each
#'  array are distributed. We overlay the density-plots on top of each array to see the
#'  distribution within the dataset.
#' @param dat input a expression data.frame (probes in row; samples in column).
#' @param annCol input an information of samples that is used to combine the sample. 
#' @param annColors input the specified colors of the sample information.
#' @param xlab character string specifying the label of X axis.
#' @param main character string specifying the graph title.    
#' @return a density-plot of microarry type of data. 
#' @seealso \code{\link{geom_density}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*6,sd=sd),100,6))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' gg_densityplot(y)
#' annCol=c("Control","Control","Control","Tumor","Tumor","Tumor" )
#' annColors <- c("Control"="black", "Tumor"="red")
#' gg_densityplot(y, annCol=annCol, annColors=annColors)
#' gg_densityplot(y, annCol=annCol)
gg_densityplot <- function(dat, annCol=NULL, annColors=NULL,
                           xlab="Expression Value (log2)",
                           main="Density-Plot"){
  df <- melt(dat, variable.name='Array')  
  if(is.null(annCol) & is.null(annColors)){
    p <- ggplot(df, aes(x=value, colour=Array)) +
      geom_density() +
      xlab(xlab) + ylab("Density") +
      labs(title=main)  
  }
  if(!is.null(annCol) & !is.null(annColors)){
    df2 <- data.frame(df, Group=rep(as.character(annCol), each=nrow(dat)))
    tmp <- annColors[annCol]
    names(tmp) <- colnames(dat)
    p <- ggplot(df2, aes(x=value, colour=Array)) +  
      geom_density() +
      scale_color_manual(values = tmp) +
      xlab(xlab) + ylab("Density") +
      labs(title=main)
  }
  if(!is.null(annCol) & is.null(annColors)){
    df2 <- data.frame(df, Group=rep(as.character(annCol), each=nrow(dat)))
    p <- ggplot(df2, aes(x=value, colour=Group)) +  
      geom_density() +
      xlab(xlab) + ylab("Density") +
      labs(title=main)
  }
  print(p)  
}


#' Generate box-plot for QC accessment on microarray type of data. 
#'
#' In processing microarry data, usually we need to see how the expression values on each
#'  array are distributed. We put box-plots side-by-side to see the distribution within
#'  the dataset.
#' @param dat input a expression data.frame (probes in row; samples in column).
#' @param annCol input an information of samples that is used to combine the sample. 
#' @param annColors input the specified colors of the sample information.
#' @param ylab character string specifying the label of y axis.
#' @param main character string specifying the graph title.    
#' @return a box-plot of microarry type of data. 
#' @seealso \code{\link{geom_boxplot}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*6,sd=sd),100,6))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' gg_boxplot(y)
#' gg_boxplot(y, x.angle=45, x.size=10)
#' gg_boxplot(y, outlier=FALSE)
#' annCol=c("Control","Control","Control","Tumor","Tumor","Tumor" )
#' annColors <- c("Control"="black", "Tumor"="red")
#' gg_boxplot(y, annCol=annCol)
#' gg_boxplot(y, annCol=annCol, annColors=annColors)
gg_boxplot <- function(dat, annCol=NULL, annColors=NULL,
                       ylab="Expression Value (log2)",
                       outlier=TRUE, x.angle=90, x.vjust=0.5, x.size=8,
                       main="Box-Plot"){
  df <- melt(dat, variable.name='Array')  
  if(is.null(annCol) & is.null(annColors)){
    if (outlier==TRUE) {
      p <- ggplot(df, aes(x=Array, y=value, fill=Array)) +
        geom_boxplot() +
        xlab("") + ylab(ylab) +
        labs(title=main)
    } else {
      p <- ggplot(df, aes(x=Array, y=value, fill=Array)) +
        geom_boxplot(outlier.shape = NA) +
        xlab("") + ylab(ylab) +
        labs(title=main) +
        scale_y_continuous(limits = quantile(df$value, c(0.1, 0.9), na.rm=TRUE), oob=scales::rescale_none)
    }
  }
  if(!is.null(annCol) & is.null(annColors)){
    df2 <- data.frame(df, Group=rep(as.character(annCol), each=nrow(dat)))
    if (outlier==TRUE) {
      p <- ggplot(df2, aes(x=Array, y=value, fill=Group)) +
        geom_boxplot() +
        #scale_fill_manual(values = annColors) +
        xlab("") + ylab(ylab) +
        labs(title=main)
    } else {
      p <- ggplot(df2, aes(x=Array, y=value, fill=Group)) +
        geom_boxplot(outlier.shape = NA) +
        #scale_fill_manual(values = annColors) +
        xlab("") + ylab(ylab) +
        labs(title=main) +
        scale_y_continuous(limits = quantile(df$value, c(0.1, 0.9), na.rm=TRUE), oob=scales::rescale_none)    
    }
  }
  if(!is.null(annCol) & !is.null(annColors)){
    df2 <- data.frame(df, Group=rep(as.character(annCol), each=nrow(dat)))
    if (outlier==TRUE) {
      p <- ggplot(df2, aes(x=Array, y=value, fill=Group)) +
        geom_boxplot() +
        scale_fill_manual(values = annColors) +
        xlab("") + ylab(ylab) +
        labs(title=main)
    } else {
      p <- ggplot(df2, aes(x=Array, y=value, fill=Group)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = annColors) +
        xlab("") + ylab(ylab) +
        labs(title=main) +
        scale_y_continuous(limits = quantile(df$value, c(0.1, 0.9), na.rm=TRUE), oob=scales::rescale_none)    
    }
  }
  print(p + theme(axis.text.x=element_text(angle=x.angle, vjust=x.vjust, size=x.size)) )
}



#' Generate principal components analysis (PCA) Plot on microarray type of data. 
#'
#' In processing microarry data, usually we need to see how the expression values on each
#'  array are separated in space.
#' @param x input a expression data.frame (probes in row; samples in column).
#' @param annCol input an information of samples that is used to distinguish the
#'  sample in colors. 
#' @param annColors input the specified colors of the sample information.
#' @param annS input an information of samples that is used to distinguish the
#'  sample in shapes. 
#' @param annShapes input the specified shapes of the sample information.
#' @param main character string specifying the graph title.    
#' @return a density-plot of microarry type of data. 
#' @seealso \code{\link{SamplePCA}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*6,sd=sd),100,6))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' annCol=c("Control","Control","Control","Tumor","Tumor","Tumor" )
#' annS=c("Batch1","Batch1","Batch2","Batch2","Batch3","Batch3" )
#' annColors <- c("Control"="black", "Tumor"="red")
#' annShapes=c("Batch1"=1,"Batch2"=4,"Batch3"=6 )
#' gg_PCA(y)
#' gg_PCA(y, annCol=annCol)
#' gg_PCA(y, annS=annS)
#' gg_PCA(y, annCol=annCol, annColors=annColors)
#' gg_PCA(y, annS=annS, annShapes=annShapes)
#' gg_PCA(y, annCol=annCol, annS=annS)
gg_PCA <- function(x, annCol=NULL, annColors=NULL,
                   annS=NULL, annShapes=NULL,
                   main="1st & 2nd PC") {
  data <- x
  spca <- ClassDiscovery::SamplePCA(data, usecor=F, center=T)
  pct1 <- round(spca@variances[1]/sum(spca@variances), digits=3)*100
  pct2 <- round(spca@variances[2]/sum(spca@variances), digits=3)*100
  xlab.text = paste("First Comp: ", as.character(pct1), "% variance", sep="")
  ylab.text = paste("Second Comp: ", as.character(pct2), "% variance", sep="")  
  
  if(is.null(annCol) & is.null(annS)){
    pc <- data.frame(spca@scores)  
    p <- ggplot(pc, aes(x=PC1, y=PC2)) +
          geom_point(size = 3)
  }
  if(!is.null(annCol) & is.null(annColors) & is.null(annS)){
    pc <- data.frame(spca@scores, Colors=annCol)  
    p <- ggplot(pc, aes(x=PC1, y=PC2)) +
          geom_point(aes(colour = Colors), size = 3)
  }
  if(!is.null(annCol) & !is.null(annColors) & is.null(annS)){
    pc <- data.frame(spca@scores, Colors=annCol)  
    p <- ggplot(pc, aes(x=PC1, y=PC2)) +
          geom_point(aes(colour = Colors), size = 3) +
          scale_colour_manual(values = annColors)    
  }
  if(is.null(annCol) & !is.null(annS) & is.null(annShapes)){
    pc <- data.frame(spca@scores, Shapes=annS)  
    p <- ggplot(pc, aes(x=PC1, y=PC2)) +
      geom_point(aes(shape = Shapes), size = 3)
  }
  if(is.null(annCol) & !is.null(annS) & !is.null(annShapes)){
    pc <- data.frame(spca@scores, Shapes=annS)  
    p <- ggplot(pc, aes(x=PC1, y=PC2)) +
      geom_point(aes(shape = Shapes), size = 3) +
      scale_shape_manual(values = annShapes)
  }
  if(!is.null(annCol) & !is.null(annS) & is.null(annColors) & is.null(annShapes)){
    pc <- data.frame(spca@scores, Colors=annCol, Shapes=annS)  
    p <- ggplot(pc, aes(x=PC1, y=PC2)) +
      geom_point(aes(colour = Colors, shape = Shapes), size = 3)
  }
  if(!is.null(annCol) & !is.null(annS) & !is.null(annColors) & is.null(annShapes)){
    pc <- data.frame(spca@scores, Colors=annCol, Shapes=annS)  
    p <- ggplot(pc, aes(x=PC1, y=PC2)) +
      geom_point(aes(colour = Colors, shape = Shapes), size = 3) +
      scale_colour_manual(values = annColors)
  }
  if(!is.null(annCol) & !is.null(annS) & !is.null(annColors) & is.null(annShapes)){
    pc <- data.frame(spca@scores, Colors=annCol, Shapes=annS)  
    p <- ggplot(pc, aes(x=PC1, y=PC2)) +
      geom_point(aes(colour = Colors, shape = Shapes), size = 3) +
      scale_colour_manual(values = annColors) +
      scale_shape_manual(values = annShapes)
  }
  final.p <- p + xlab(xlab.text) + ylab(ylab.text) + labs(title=main) +
    geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") +
    geom_vline(aes(xintercept=0), colour="#990000", linetype="dashed") 
  print(final.p) 
}  



#' Generate correlation-plot between X and Y.
#'
#' Generate correlation-plot between X and Y, returning results for both Pearson and
#'  Spearman correlation. 
#' @param x input numeric vector to specify the coordinate on X-axis.
#' @param y input numeric vector to specify the coordinate on Y-axis.
#' @param xlab character string specifying the label of X axis.
#' @param ylab character string specifying the label of y axis.
#' @param dot.size input the size of the dot. Default is 2.     
#' @return a correlation-plot of X and Y. 
#' @seealso \code{\link{rcorr_XY}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*6,sd=sd),100,6))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' gg_corplot_XY(x=y[,1], y=y[,2])
#' gg_corplot_XY(x=y[,1], y=y[,2], xlab="Random X", ylab="Random Y")
gg_corplot_XY <- function(x, y, xlab=NULL, ylab=NULL, dot.size=2){
  if(is.null(xlab)) xlab <- deparse(substitute(x))
  if(is.null(ylab)) ylab <- deparse(substitute(y))
  ### combine data ###
  dat <- data.frame(x=x, y=y)
  
  ### calculate both pearson and spearman correlation
  cor_p <- rcorr_XY(dat$x, y=dat$y, method='pearson')
  cor_s <- rcorr_XY(dat$x, y=dat$y, method='spearman')
  
  gtitle <- sprintf('Correlation \n Pearson - r.coef: %s, p-val: %s\n Spearman - rho.coef: %s, p-val: %s\n',paste(round(cor_p$r.coef, 3), collapse=', '),paste(round(cor_p$pearson.pval, 5), collapse=', '),paste(round(cor_s$rho.coef, 3), collapse=', '),paste(round(cor_s$spearman.pval, 5) ,collapse=', '))
  toplot <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(size=dot.size) + 
    geom_smooth(method = "lm", se = TRUE) +
    labs(x=xlab, y=ylab,
         title=gtitle) +
    theme_bw()
  return(toplot)
} 



#' Generate concordance-plot between X and Y.
#'
#' Generate concordance-plot between X and Y, returning concordance coefficient from 
#'  'Preprocess' package in R.
#' @param x input numeric vector to specify the coordinate on X-axis.
#' @param y input numeric vector to specify the coordinate on Y-axis.
#' @param xlab character string specifying the label of X axis.
#' @param ylab character string specifying the label of y axis.
#' @param dot.size input the size of the dot. Default is 2.    
#' @return a concordance-plot of X and Y. 
#' @seealso \code{\link{f.cord}} which this function wraps.
#' @export
#' @examples
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*6,sd=sd),100,6))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' gg_concord_XY(x=y[,1], y=y[,2])
#' gg_concord_XY(x=y[,1], y=y[,2], xlab="Random X", ylab="Random Y")
gg_concord_XY <- function(x, y, xlab=NULL, ylab=NULL, dot.size=2){
  if(is.null(xlab)) xlab <- deparse(substitute(x))
  if(is.null(ylab)) ylab <- deparse(substitute(y))
  ### combine data ###
  dat <- data.frame(x=x, y=y)
  ### calculate both pearson and spearman correlation
  cord_p <- PreProcess::f.cord(dat$x, dat$y)
  
  gtitle <- sprintf('Concordance \n coefficient: %s',paste(round(cord_p, 3), collapse=', '))
  toplot <- ggplot(dat, aes(x = x, y = y)) +
    geom_abline(colour="red", size=1) +
    geom_point(size=dot.size) + 
    labs(x=xlab, y=ylab,
         title=gtitle) +
    theme_bw()
  return(toplot)
} 





#' Generate bar-plot for the frequency of elements in data.frame. 
#'
#' In processing microarry data, usually we need to see how the expression values on each
#'  array are distributed. We put box-plots side-by-side to see the distribution within
#'  the dataset.
#' @param dat input a data.frame (probes in row; samples in column).
#' @param percentage Default is FALSE. It will show the counts of the element
#'  when it is FALSE.
#' @param annColors input the specified colors of the element.
#' @param label  Default is FALSE. When it is TRUE, it will show the 
#'  count/percentage in the figure.
#' @param ylim input the limits in y axis.
#' @param x_order input the order in x axis. 
#' @param xlab character string specifying the label of x axis.
#' @param ylab character string specifying the label of y axis.
#' @param main character string specifying the graph title.
#' @param legend.title input the title of the legend.    
#' @return a bar-plot of a data frame. 
#' @seealso \code{\link{geom_bar}} which this function wraps.
#' @export
#' @examples
#' y <- data.frame(matrix(0, 9, 10))
#' rownames(y) <- paste("Probe",1:9, sep="_")
#' colnames(y) <- c(paste("A", 1:10, sep=""))
#' dat <- data.frame(apply(y, c(1,2), function(x) sample(c("Gain","Loss", NA),1)))
#' cols <- c("Gain"="dodgerblue2", "Loss"="indianred2")
#' gg_order <- rownames(dat)[order(apply(apply(dat, 1, table),2, sum), decreasing = TRUE)]
#' #gg_order <- rownames(dat)[order(do.call("rbind", lapply(apply(dat, 1, table), sum)), decreasing = TRUE)]
#' gg_frequency(dat)
#' gg_frequency(dat, percentage=TRUE)
#' gg_frequency(dat, percentage=TRUE, annColors=cols)
#' gg_frequency(dat, percentage=TRUE, annColors=cols, label=TRUE)
#' gg_frequency(dat, percentage=TRUE, annColors=cols, label=TRUE, ylim=c(0, 0.9))
#' gg_frequency(dat, percentage=TRUE, annColors=cols, label=TRUE, ylim=c(0, 0.9),
#'              x_order=gg_order)
#' gg_frequency(dat, percentage=TRUE, annColors=cols, label=TRUE, ylim=c(0, 0.9),x_order=gg_order,
#'              ylab="Percentage of Alteration", x.angle=0, legend.title="Alteration")
#' gg_frequency(dat, annColors=cols)
#' gg_frequency(dat, annColors=cols, label=TRUE)
#' gg_frequency(dat, annColors=cols, label=TRUE, ylim=c(0, 30))
#' gg_frequency(dat, annColors=cols, label=TRUE, ylim=c(0, 30),
#'              gg_order=gg_order)
#' gg_frequency(dat, annColors=cols, label=TRUE, ylim=c(0, 30), x_order=gg_order,
#'              ylab="Percentage of Alteration", x.angle=0, legend.title="Alteration")
gg_frequency <- function(dat, percentage=FALSE, annColors=NULL, label=FALSE,
                         ylim=NULL, x_order=NULL,
                         xlab="",
                         ylab="Count/Percentage",
                         legend.title=NULL,
                         x.angle=90, x.vjust=0.5, x.size=8,
                         main=""){
  tmp <- na.omit(melt(as.matrix(dat)))
  tmp2 <- table(tmp$Var1, tmp$value)
  count <- melt(tmp2)
  freq <- count
  freq$value <- count$value/ncol(dat)
  
  ### remove if value==0
  if(length(which(count$value == 0))!=0){
    count <- count[-which(count$value == 0),]
    freq <- freq[-which(freq$value == 0),]
  }
  
  if (is.null(x_order)==FALSE){
    count$Var1 <- factor(count$Var1, levels=x_order)
    freq$Var1 <- factor(freq$Var1, levels=x_order)
  }
  if (percentage==TRUE){
    p <- ggplot(freq, aes(x = Var1, y = value, fill = Var2,
                          label=paste(round(value*100), "%", sep=""))) + 
            geom_bar(stat = "identity") 
  }
  if (percentage==FALSE){
    p <- ggplot(count, aes(x = Var1, y = value, fill = Var2,
                          label=value)) + 
            geom_bar(stat = "identity")
  }
  if (is.null(annColors)==FALSE){
    p <- p + scale_fill_manual(values=cols)
    #scale_fill_manual(values=cols, breaks=c("Gain", "Loss")) +
  }
  if (label==TRUE){
    p <- p + geom_text(position="stack", vjust=0)
    #geom_text(position="stack", aes(ymax=1),vjust=0)
  }
  if (is.null(ylim)==TRUE){
    if (percentage==TRUE){
      p <- p + scale_y_continuous(labels = percent)
    }
    if (percentage==FALSE){
      p <- p + scale_y_continuous(limits=ylim, oob=scales::rescale_none)
    } 
  }
  if (is.null(ylim)==FALSE){
    if (percentage==TRUE){
      p <- p + scale_y_continuous(labels = percent, limits=ylim, oob=scales::rescale_none)
    }
    if (percentage==FALSE){
      p <- p + scale_y_continuous(limits=ylim, oob=scales::rescale_none)
    } 
  }
  return(p + xlab(xlab) + ylab(ylab) + labs(title=main) +
         theme(axis.text.x=element_text(angle=x.angle, vjust=x.vjust, size=x.size)) +
         guides(fill = guide_legend(title = legend.title)) )
}
  
  

#' Make a Venn Diagram with count table that each component corresponding to a separate circle in VD.
#'
#' This function takes a list and creates a publication-quality Venn Diagram. 
#' @param x A list of vectors (e.g., integers, chars), with each component corresponding
#'  to a separate circle in the Venn diagram
#' @param filename Default is set to NULL which returns the grid object itself.
#' @param height Integer giving the height of the output figure in units.
#' @param width  Integer giving the width of the output figure in units.
#' @param resolution Resolution of the final figure in DPI.
#' @param units Size-units to use for the final figure.
#' @param compression What compression algorithm should be applied to the final tiff.
#' @param na Missing value handling method: "none", "stop", "remove".
#' @param main Character giving the main title of the diagram.
#' @param sub Character giving the subtitle of the diagram. 
#' @param main.pos Vector of length 2 indicating (x,y) of the main title.
#' @param main.fontface Character giving the fontface (font style) of the main title.
#' @param main.fontfamily Character giving the fontfamily (font type) of the main title.
#' @param main.col Character giving the colour of the main title.
#' @param main.cex Number giving the cex (font size) of the main title.
#' @param main.just Vector of length 2 indicating horizontal and vertical justification of the main title.
#' @param sub.pos Vector of length 2 indicating (x,y) of the subtitle.
#' @param sub.fontface Character giving the fontface (font style) of the subtitle.
#' @param sub.fontfamily Character giving the fontfamily (font type) of the subtitle.
#' @param sub.col Character Colour of the subtitle.
#' @param sub.cex Number giving the cex (font size) of the subtitle.
#' @param sub.just Vector of length 2 indicating horizontal and vertical justification of the subtitle.
#' @param category.names Allow specification of category names using plotmath syntax.
#' @param force.unique Logical specifying whether to use only unique elements in each item of the input
#'  list or use all elements. Defaults to FALSE
#' @param ... A series of graphical parameters tweaking the plot. See below for details.  
#' @return A Venn Diagram with a count table that each component corresponding to a separate circle in VD.
#' @seealso \code{\link{venn.diagram}} which this function wraps.
#' @export
#' @examples
#' mylist <- list(A = 1:150, B = 121:170)
#'gg_venn.diagram(x=mylist,
#'                col = "transparent",
#'                fill=rainbow(2), alpha=0.5,
#'                cex = 1.5, cat.fontface = 4, lty =2, 
#'                fontfamily =3,
#'                margin = 0.1, main="mylist",
#'                countFile="vennDiag_mylist.csv")
#'mylist <- list(A = 1:150, B = 121:170, C = 101:200)
#'gg_venn.diagram(x=mylist,
#'                col = "transparent",
#'                fill=rainbow(3), alpha=0.5,
#'                cex = 1.5, cat.fontface = 4, lty =2, 
#'                fontfamily =3,
#'                margin = 0.1, main="mylist",
#'                countFile="vennDiag_mylist.csv")
#'mylist <- list(
#'  I = c(1:60, 61:105, 106:140, 141:160, 166:175, 176:180, 181:205, 206:220),
#'  IV = c(531:605, 476:530, 336:375, 376:405, 181:205, 206:220, 166:175, 176:180),
#'  II = c(61:105, 106:140, 181:205, 206:220, 221:285, 286:335, 336:375, 376:405),
#'  III = c(406:475, 286:335, 106:140, 141:160, 166:175, 181:205, 336:375, 476:530)
#')
#'gg_venn.diagram(x=mylist,
#'                col = "transparent",
#'                fill=rainbow(4), alpha=0.5,
#'                cex = 1.5, cat.fontface = 4, lty =2, 
#'                fontfamily =3,
#'                margin = 0.1, main="mylist",
#'                countFile="vennDiag_mylist.csv")
gg_venn.diagram <- function(x, countFile="vennDiag_counts.csv", ...){
  if (is.list(x)==FALSE) {
    stop("x must be a list of vectors.\n")
  }
  ### venn diagram:
  v1 <- VennDiagram::venn.diagram(x, filename =NULL, ...)
  grid.newpage()
  grid.draw(v1)
  
  ### counts:
  universe.list <- sort(Reduce(union, x))
  Counts <- matrix(0, nrow=length(universe.list), ncol=length(x))
  colnames(Counts) <- names(x)
  rownames(Counts) <- universe.list
  
  for (i1 in 1:nrow(Counts)){
    for (i2 in 1:ncol(Counts)){
      Counts[i1, i2] <- universe.list[i1] %in% x[[i2]]
    }
  }
  write.csv(Counts, file=countFile)
}



#' Generate Quantile-Quantile Plots.
#'
#' Generate a normal QQ plot of the values in vec using qqnorm function and adds a line
#'  to a "theoretical", by default normal, quantile-quantile plot which passes through
#'  the probs quantiles, by default the first and third quartiles.
#' @param vec input numeric vector.
#' @param dot.size input the size of the dot. Default is 3.
#' @param xlab character string specifying the label of X axis.
#' @param ylab character string specifying the label of y axis.
#' @param main Character giving the main title of the diagram.     
#' @return a correlation-plot of X and Y. 
#' @seealso \code{\link{qqnorm}} and \code{\link{qqline}} which this function wraps.
#' @export
#' @examples
#' vec <- rt(200, df = 5)
#' gg_qqplot(vec=vec)
#' gg_qqplot(vec=vec, size=2)
gg_qqplot <- function(vec, dot.size=3, xlab="Theoretical Quantiles",
                      ylab="Sample Quantiles",
                      main="Normal Q-Q Plot"){
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  xq <- qqnorm(vec, plot.it = FALSE)$x
  yq <- qqnorm(vec, plot.it = FALSE)$y
  
  d <- data.frame(Theoretical=xq, Sample=yq)
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
                                   "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  if (var(yq)!=0){
    p <- ggplot(d, aes(x=Theoretical, y=Sample, color=Sample))
    return(ggplot(d, aes(x=Theoretical, y=Sample, color=Sample)) + geom_point(size = dot.size) +
      geom_abline(slope = slope, intercept = int) +
      xlab(xlab) + ylab(ylab) + labs(title=main) +
      scale_colour_gradientn(colours = jet.colors(7)))
  }
  if (var(yq)==0){
    p <- ggplot(d, aes(x=Theoretical, y=Sample, color="red"))
  }
  return(p + geom_point(size = dot.size) +
           geom_abline(slope = slope, intercept = int) +
           xlab(xlab) + ylab(ylab) + labs(title=main))
}




