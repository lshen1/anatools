% Generated by roxygen2 (4.0.2.9000): do not edit by hand
% Please edit documentation in R/utility_function.R
\name{get_robj}
\alias{get_robj}
\title{Get one of the objects from a R file/RData.}
\usage{
get_robj(x, RData)
}
\arguments{
\item{x}{input the name of object that users specify.}

\item{RData}{input a R file/RData containing objects.}
}
\value{
an object that user specified from a R file/RData.
}
\description{
This function allows you to get one of the objects from
 a R file/RData.
}
\examples{
x=1:10
y=letters[1:10]
save(x, y, file="xy.RData")
dat.dir <- file.path("./", "xy.RData")
show_robj(dat.dir)
rm(x, y)
x <- get_robj(x="x", dat.dir)
x
}
\seealso{
\code{\link{get}} which this function wraps.
}

