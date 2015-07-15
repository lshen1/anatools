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

