% Generated by roxygen2 (4.0.2.9000): do not edit by hand
% Please edit documentation in R/autoplot.R
\name{autoplot.tableAndPlot}
\alias{autoplot.tableAndPlot}
\title{Arrange and plot a survival plot, it's legend and a table.}
\usage{
\method{autoplot}{tableAndPlot}(object, ..., hideTabLeg = TRUE,
  plotHeight = 0.8, tabHeight = 0.2)
}
\arguments{
\item{object}{An object of class \code{tableAndPlot} as returned by
\code{autoplot.survfit}}

\item{...}{Additional arguments (not implemented)}

\item{hideTabLeg}{Suppress table legend. If \code{supTabLeg = FALSE}
then a legend will also be shown for the table}

\item{plotHeight}{Height of plot}

\item{tabHeight}{Height of table}
}
\value{
A graph, as plotted by \code{gridExtra::grid.arrange}
}
\description{
Uses \code{gridExtra::gridArrange}
to arrange a plot, it's legend and a table.
}
\details{
Put Kaplan-Meier Plot and table together using `survMisc` package in R

Put Kaplan-Meier Plot and table together using `survMisc` package in R.

Arguments to \code{plotHeigth} and \code{tabHeight} are
best specified as fractions adding to one e.g. \eqn{0.85 + 0.15 =1}.
\cr \cr
Other \code{ggplot2::ggplot2} objects may be plotted using this
method. They need to be stored in a list of two (with element \code{table}
if \code{hideTabLeg} is to be used). The class of this list should be
modified with
\cr \cr
\code{class(list1) <- c("tableAndPlot", "list")}
}
\examples{
data(colon)
fit <- survfit(Surv(time,status)~rx, data=colon)
autoplot(autoplot(fit, timeby=500, xlab="Time (Day)"))
mycol <- c("rx=Lev"="blue", "rx=Lev+5FU"="red", "rx=Obs"="black")
autoplot(autoplot(fit, timeby=500, xlab="Time (Day)", assignCol=mycol))
}
\author{
Chris Dardis. Based on existing work by
R. Saccilotto, Abhijit Dasgupta, Gil Tomas and Mark Cowley.
}

