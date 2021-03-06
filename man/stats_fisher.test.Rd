% Generated by roxygen2 (4.0.2.9000): do not edit by hand
% Please edit documentation in R/modelTesting.R
\name{stats_fisher.test}
\alias{stats_fisher.test}
\title{Perform Fisher's exact test between each row in a matrix and a vector using `stats` package in R}
\usage{
stats_fisher.test(data, group, rowlevel, ...)
}
\arguments{
\item{data}{a matrix with at least 2 outcomes (0/1; Mut/WT, etc).}

\item{group}{a vector of factor which will be concatenated to the row of 'data'.}

\item{rowlevel}{specifies the order of factor use in a contingency table.}

\item{...}{other parameters in fisher.test function.}
}
\value{
the contingency table and p-value.
}
\description{
Performs Fisher's exact test for testing the null of independence of rows and columns
 in a contingency table with fixed marginals.
}
\examples{
m0 <- matrix(0, 5, 30)
dat <- apply(m0, c(1,2), function(x) sample(c(0,1),1))
rownames(dat) <- paste("R", 1:5, sep="")
colnames(dat) <- paste("C", 1:30, sep="")
si.gp <- factor(c(rep("E", 15), rep("M", 15)), levels=c("E", "M"))
stats_fisher.test(dat, si.gp, rowlevel=c("0", "1"))
}
\seealso{
\code{\link{fisher.test}} which this function wraps.
}

