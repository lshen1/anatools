#' Show a list of names of objects from a R file/RData.  
#'
#' This function allows you to list the names of objects from
#'  a R file/RData
#'  
#' @param RData input a R file/RData containing objects.
#' 
#' @return a vector of character strings giving the names of
#'  the objects in a R file/RData.
#'  
#' @seealso \code{\link{objects}} which this function wraps.
#' @export
#' @examples
#' x=1:10
#' y=letters[1:10]
#' save(x, y, file="xy.RData")
#' dat.dir <- file.path("./", "xy.RData")
#' show_robj(dat.dir)
show_robj <- function(RData){
  tempEnv <- new.env()
	load(RData, envir=tempEnv)
	return(objects(tempEnv))
}


#' Get one of the objects from a R file/RData.  
#'
#' This function allows you to get one of the objects from
#'  a R file/RData.
#'  
#' @param x input the name of object that users specify.    
#' @param RData input a R file/RData containing objects.
#' 
#' @return an object that user specified from a R file/RData.
#' 
#' @seealso \code{\link{get}} which this function wraps.
#' @export
#' @examples
#' x=1:10
#' y=letters[1:10]
#' save(x, y, file="xy.RData")
#' dat.dir <- file.path("./", "xy.RData")
#' show_robj(dat.dir)
#' rm(x, y)
#' x <- get_robj(x="x", dat.dir)
#' x
get_robj <- function(x, RData){
  tempEnv <- new.env()
  load(RData, tempEnv)
  res <- get(x, envir=tempEnv)
  return(res)
}

#' Is an Object Defined?  
#'
#' Look for an R object of the given name.
#' 
#' @param object input a variable name (doesn't need to be a character string).
#'     
#' @return Logical, true if and only if an object of the correct name
#'  and mode is found.
#'  
#' @seealso \code{\link{exists}} which this function wraps.
#' @export
#' @examples
#' x=1:10
#' exists_robj(x)
#' exists("x")  #using exists, the object must be a character string
exists_robj <- function(object)
{
  return(exists(as.character(substitute(object))))
}



#' Aggregate the microarry data based on the maximum variance probe of gene probesets. 
#'
#' In microarray analysis, some genes are associated with multiple probes. In practice,
#'  sometime we need to summarize the expression in gene level (each gene has only one
#'  value of a given sample). This function would aggregate the microarry data based on
#'  the maximum variance probe of gene probesets.
#'    
#' @param dat input a expression data.frame (probes in row; samples in column).  
#' @param annot input the annotation of the expression data (must have identical rownames).
#' @param symcol default is set to "Symbol" which defines the column with multiple probes
#'  information from the annotation.
#'      
#' @return a list contains the list of "data" and "annot" with maximum variance probesets.
#'  
#' @seealso \code{\link{var}} which this function wraps.
#' @export
#' @examples
#' # Simulate gene expression data for 100 probes and 6 microarrays
#' # Microarray are in two groups
#' # First 50 probes are differentially expressed in second group
#' # Std deviations vary between genes with prior df=4
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*6,sd=sd),100,6))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' set.seed(123)
#' annot.y <- data.frame(ID=rownames(y),
#'  Symbol=paste("Gene", sort(sample(seq(1, 75), 100, replace=TRUE)), sep="_"))
#' rownames(annot.y) <- rownames(y)
#' mvDat <- aggregate_maxvariance(y, annot.y, symcol="Symbol")
#' names(mvDat)                     
aggregate_maxvariance <- function(dat, annot, symcol="Symbol"){
  if (identical(rownames(dat), rownames(annot))==FALSE){
    stop("Data and Annotation must have identical rownames!")
  } else {
    usedat <- data.frame(var=apply(dat, 1, var, na.rm=TRUE), Symbol=annot[,symcol],
                         stringsAsFactors=FALSE)
    dat.ls <- split(usedat, usedat$Symbol)
    ann.ls <- split(annot, annot[,symcol])
    ind.ls <- lapply(unname(dat.ls), function(x) x[which.max(x$var),])
    ind <- rownames(do.call("rbind", ind.ls))
    res <- list(data=dat[ind,], annot=annot[ind,])
    return(res)
  }
}


#' Aggregate the microarry data based on the mean of gene probesets. 
#'
#' In microarray analysis, some genes are associated with multiple probes. In practice,
#'  sometime we need to summarize the expression in gene level (each gene has only one
#'  value of a given sample). This function would aggregate the microarry data based on
#'  the mean of gene probesets.
#'    
#' @param dat input a expression data.frame (probes in row; samples in column).  
#' @param annot input the annotation of the expression data (must have identical rownames).
#' @param symcol default is set to "Symbol" which defines the column with multiple probes
#'  information from the annotation.
#'      
#' @return a list contains the aggregated "data" and "annot".
#'  
#' @seealso \code{\link{aggregate}} which this function wraps.
#' @export
#' @examples
#' # Simulate gene expression data for 100 probes and 6 microarrays
#' # Microarray are in two groups
#' # First 50 probes are differentially expressed in second group
#' # Std deviations vary between genes with prior df=4
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- data.frame(matrix(rnorm(100*6,sd=sd),100,6))
#' rownames(y) <- paste("Probe",1:100, sep="_")
#' colnames(y) <- c(paste("A", 1:3, sep=""), paste("B", 1:3, sep=""))
#' y[1:50,4:6] <- y[1:50,4:6] + 2
#' y[sample(seq(1,100), 10), 1] <- NA
#' y[sample(seq(1,100), 10), 3] <- NA
#' y[sample(seq(1,100), 10), 5] <- NA
#' set.seed(123)
#' annot.y <- data.frame(ID=rownames(y),
#'  Symbol=paste("Gene", sort(sample(seq(1, 75), 100, replace=TRUE)), sep="_"))
#' rownames(annot.y) <- rownames(y)
#' mnDat <- aggregate_mean(y, annot.y, symcol="Symbol")
#' names(mnDat)                     
aggregate_mean <- function(dat, annot, symcol="Symbol"){
  if (identical(rownames(dat), rownames(annot))==FALSE){
    stop("Data and Annotation must have identical rownames!")
  } else {
    dat.avg <- aggregate(dat, list(Symbol=annot[,symcol]), mean, na.rm=TRUE)
    rownames(dat.avg) <- dat.avg[, 1]
    dat.avg <- dat.avg[, -1]
    ind <- rownames(dat.avg)
    ann.avg <- annot[match(ind, annot[,symcol]),]
    rownames(ann.avg) <- ann.avg[,symcol]
    res <- list(data=dat.avg, annot=ann.avg)
    return(res)
  }
}



#' Create directory/folder after checking directory/folder existence using a direct path 
#'
#' This function would check the directory existence. If the folder doesn't exist,
#'  it creates the folder with permission mode to `777`.
#'  
#' @param fname input the folder user would like to create. 
#'       
#' @return a folder with user defined name if the folder doesn't exist.
#' 
#' @seealso \code{\link{dir.create}} which this function wraps. Also see \code{\link{create_folders}}
#'  for using a character string.
#' @export
#' @examples
#' create_folder("./test")
#' list.dirs("./", full.names=FALSE, recursive = FALSE)
#' fName <- "./test2"
#' create_folders(fName)
#' list.dirs("./", full.names=FALSE, recursive = FALSE)
create_folder <- function(fname){
  if (!file.exists(as.character(substitute(fname)))){
    dir.create(as.character(substitute(fname))) 
  }
}

#' Create directory/folder after checking directory/folder existence using a character string 
#'
#' This function would check the directory existence. If the folder doesn't exist,
#'  it creates the folder with permission mode to `777`.
#'  
#' @param fname input the folder user would like to create.
#'        
#' @return a folder with user defined name if the folder doesn't exist.
#' 
#' @seealso \code{\link{dir.create}} which this function wraps. Also see \code{\link{create_folder}}
#'  for using direct path.
#' @export
#' @examples
#' create_folder("./test")
#' list.dirs("./", full.names=FALSE, recursive = FALSE)
#' fName <- "./test2"
#' create_folders(fName)
#' list.dirs("./", full.names=FALSE, recursive = FALSE)
create_folders <- function(fname){
  if (!file.exists(fname)){
    dir.create(fname) 
  }
}


#' Check availability of a gene list in the annotation of the microarry data. 
#'
#' In microarray data analysis, sometimes we only use a set of genes that PI selected for
#'  further analysis. We would like to first see the availability of the gene list in
#'  the annotation of the microarray data.
#'  
#' @param Annot input the vector of characters in annotation used to match the gene list.
#' @param GeneList input a gene list.
#'        
#' @return a summary of probes of the gene list in annotation.
#' 
#' @export
#' @examples
#' x <- c(sort(sample(3:23, 100, replace=TRUE)), NA)
#' y <- c(sort(sample(1:20, 9)), NA)
#' avail_SIGinANN(x, y)
avail_SIGinANN <- function(Annot, GeneList){
  ### common features in both sets ###
  mSig <- intersect(Annot, GeneList)
  cat("------------------------------------------------------------------", "\n")
  cat("Common Features In Both Annotation and Genelist:", "\n")
  print(length(mSig))
  
  ### The number of probesets can be match for GeneList ###
  cat("------------------------------------------------------------------", "\n")
  cat("The number of probesets can be found in Annotation for Genelist:", "\n")
  print(sort(table(as.character(Annot[which(Annot %in% mSig)])), decreasing = TRUE))  
  
  ### Symbols cannot be matched ###
  cat("------------------------------------------------------------------", "\n")
  cat("Which feature in the GeneList cannot be found in Annotation:", "\n")
  print(GeneList[which(!GeneList %in% mSig)])
}



