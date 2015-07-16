#' Generate Kaplan-Meier Plot and table using `survMisc` package in R 
#'
#' Generate Kaplan-Meier Plot and table using `survMisc` package in R.
#' 
#' @name autoplot.survfit
#' @export autoplot.survfit
#' @import data.table
#' @method autoplot survfit
#' @description Uses \code{ggplot2::ggplot2} to generate survival curves
#'  (Kaplan-Meier plot) and a table showing no. of events per time period
#' @title Generate a ggplot2::ggplot for \code{survfit} object
#' 
#' @param object An object of class \code{survfit}
#' @param ... Additional arguments (not implemented)
#' @param timeby the survival time interval 
#' @param xlab Label for x axis on survival plot
#' @param ylab Label for y axis on survival plot
#' @param title Title for survival plot
#' @param titTextSize Title size for survival plot
#' @param axisTitSize Title size for axes
#' @param axisLabSize Title size for label axes
#' @param survSize Survival line size
#' @param assignCol default is NULL
#' @param type If \code{type="single"} (the default), plots single lines.
#' \cr \cr
#' If \code{type="CI"} will add lines indicating confidence intervals
#' (taken from \code{upper} and \code{lower} values of \code{survfit} object.
#' Higher values of \code{alpha} (tansparency) are recommended
#' for this, e.g. \code{alpha=0.8}.
#' \cr \cr
#' If \code{type="fill"} will add filled rectangles from the survival lines to
#' the confidence intervals above.
#' @param palette Options are taken from
#' \href{http://colorbrewer2.org/}{color_brewer}.
#' \cr \cr
#' \code{palette="Dark2"} (the default) is recommended for
#' ##' \code{single} or \code{CI} plots.
#' \cr \cr
#' \code{palette="Set2"} is recommended for \code{fill} plots.
#' @param jitter If \code{jitter="noEvents"}, adds some random, positive noise
#' to survival lines with no events (i.e. all observations censored).
#' This will bring them just above 1 on the y-axis, making them easier
#' to see separately.
#' \cr \cr
#' If \code{jitter="all"} add some vertical noise to all survival lines
#' @param legend If \code{legend=FALSE}, no legends will be produced
#' for the plot or table
#' @param legLabs These can be used to replace the names
#' of the strata from the fit. Should be given in the same
#' order as those strata
#' @param legTitle Title for legend
#' @param legTextSize Title size for legend
#' @param legSize Legend (key) width and height
#' @param alpha Alpha, transparency of lines indicating confidence intervals
#' or filled rectangles. Should be in range \eqn{0-1}.
#' Larger values e.g. \code{alpha=0.7} are recommended for confidence
#' intervals
#' @param censShape Shape of marks to indicate censored onservations.
#' \cr Default is 3 which gives vertical ticks
#' \cr Use 10 for circular marks
#' @param censSize Size of marks to indicate censored onservations
#' @param leg.x input the x.axis position of the legend
#' @param leg.y input the y.axis position of the legend
#' @param CIline Confidence interval line type
#' @param fillLineSize Line size surrouding filled boxes
#' @param pval If \code{pval=TRUE}, adds \eqn{p} value from
#' log-rank test to plot
#' @param sigP No. of significant digits to display in \eqn{p} value.
#' Typically \eqn{1-3}
#' @param pX Location of \eqn{p} value on x axis. Should be in range of
#' \eqn{0 - 1}, where value is to be placed relative to the maximum observed
#' time. E.g. \code{pX = 0.5} will place it half-way along x-axis
#' @param pY Location of \eqn{p} value on y axis. Should be in range of
#' \eqn{0 - 1}, as above
#' @param tabTime If \code{tabTime="major"} this will use the major
#' x-axis (time) marks from the survival plot.
#' \cr
#' If \code{tabTime="minor"}, minor axis marks are used instead
#' @param tabTitle Table title
#' @param tabTitTextSize Table title text size
#' @param tabLegTextSize Table legend text size
#' @param nRiskSize No. at risk - text size
#' @return A \code{list} of \code{ggplot2::ggplot} objects, with elments \code{plot},
#' the survial plot and \code{table} the table of events per time.
#' This \code{list} has the additional \code{class} of
#' \code{tableAndPlot}, allowing methods from \code{autoplot.tableAndPlot}.
#' Further modifications may be made to the objects in the list if desired.
#' @author Chris Dardis. Based on existing work by
#' R. Saccilotto, Abhijit Dasgupta, Gil Tomas and Mark Cowley.
#' @examples
#' data(colon)
#' fit <- survfit(Surv(time,status)~rx, data=colon)
#' autoplot(fit)
#' autoplot(fit, type="CI", pval=TRUE, pX=0.3,
#'  title="Time to infection following catheter placement \n
#'  by type of catheter, for dialysis patients")$plot
#' fit <- survfit(Surv(time,status)~1, data=colon)
#' autoplot(fit, legLabs="")$plot
#' autoplot(fit, legend=FALSE, type="fill")$plot
autoplot.survfit <- function(object, ...,
                             timeby=500,
                             xlab="Time",
                             ylab="Survival Probability",
                             title="Kaplan-Meier Plot",
                             titTextSize=15,
                             axisTitSize=12,
                             axisLabSize=10,
                             survSize=1.0,
                             assignCol=NULL,
                             type=c("single", "CI", "fill"),
                             palette=c("Dark2","Set2","Accent","Paired",
                                       "Pastel1","Set1","Set3"),
                             jitter=c("none", "noEvents", "all"),
                             censShape=3,
                             censSize=3,
                             leg.x=0.8,
                             leg.y=0.85,
                             legend=FALSE,
                             legLabs=NULL,
                             legTitle="Strata",
                             legTextSize=10,
                             legSize=1,
                             alpha=0.05,
                             CIline=10,
                             fillLineSize=0.05,
                             pval=TRUE,
                             sigP=3,
                             pX=0.1,
                             pY=0.1,
                             tabTime=c("major", "minor"),
                             tabTitle="Number at risk by time",
                             tabTitTextSize=15,
                             tabLegTextSize=5,
                             nRiskSize=3.5){
  stopifnot(inherits(object, "survfit"))
  if(!is.null(legLabs) &! length(object$strata)==0) stopifnot(
    length(legLabs)==length(object$strata))
  ### generate data to plot
  ### declare variables (for R CMD check)
  ### st1 is vector for strata identification
  surv <- n.risk <- n.censor <- n.event <- upper <- lower <- NULL
  .SD <- st1 <- stNames <- st <- s1 <- minT <- l <- maxT <- u <- NULL
  ### change names for strata to legLabs if required
  if(is.null(legLabs)){
    stNames <- names(object$strata)
  } else {
    stNames <- legLabs
  }
  ### define breaks
  times <- seq(0, max(object$time), by = timeby)
  ### if only one strata (intercept only model)
  if (is.null(object$strata)) {
    if(is.null(legLabs)) {
      st1 <- as.factor(rep(1, length(object$time)))
    } else {
      stopifnot(length(legLabs)==1)
      st1 <- as.factor(rep(legLabs, length(object$time)))
    }
  } else {
    ### add vector for one strata according to number of rows of strata
    st1 <- unlist(sapply( 1:length(object$strata),
                          function (i) rep(stNames[i], object$strata[i]) ))
  }
  ### create data.table with data from survfit
  ### add column for strata
  ### (using data.table here as avoids duplication when adding rows later)
  ### also rename strata as 'st' to avoid calling survival::function
  dt1 <- data.table::data.table(time=object$time,
                                n.risk=object$n.risk,
                                n.event=object$n.event,
                                n.censor = object$n.censor,
                                surv=object$surv,
                                upper=object$upper,
                                lower=object$lower,
                                st=as.factor(st1))
  ### make two rows for each stratum
  ### for time=0 to time=time of first event
  dt2 <- data.table::rbindlist(list(dt1[, .SD[1, ], by=st],
                        dt1[, .SD[1, ], by=st]))
  ### set n.event and n.censored to zero
  dt2[, c("n.event", "n.censor") := list(0), by=st]
  ### set surv, upper and lower to one
  dt2[, c("surv", "upper", "lower") := list(1), by=st]
  ### set first time to zero
  dt2[seq(length(unique(dt2$st))), "time" := (0L) ]
  ### reorder to allow binding
  data.table::setcolorder(dt2, names(dt1))
  dt1 <- data.table::rbindlist(list(dt2, dt1))
  ### jitter
  jitter <- match.arg(jitter)
  ### for groups with no events add random no.to survival (by strata)
  if (jitter=="noEvents") {
    ### add column to indicate no. events by group
    dt1[, s1 := sum(n.event), by=list(st)]
    dt1[s1==0, surv := surv+(runif(1, 0.01, 0.05)), by=st]
  }
  if(jitter=="all"){
    ### for groups with no events add random no.to survival (by strata)
    dt1[, surv := surv+(runif(1, 0.01, 0.05)), by=st]
  }
  ###
  dt1 <- dt1[order(st)]
  ### plot single lines only
  g1 <- ggplot2::ggplot(data=dt1, ggplot2::aes(group=st, colour=st)) +
    ggplot2::geom_step(ggplot2::aes(x=time, y=surv), direction="hv", size=survSize) + ggplot2::theme_bw()
  ###
  type <- match.arg(type)
  if (type=="CI"){
    g1 <- g1 +
      ggplot2::geom_step(ggplot2::aes(x=time, y=upper),
                direction="hv", linetype=CIline, alpha=alpha) +
      ggplot2::geom_step(ggplot2::aes(x=time, y=lower),
                direction="hv", linetype=CIline, alpha=alpha)
  }
  if (type=="fill"){
    ### copy dt1 to work allow further work
    dt2 <- dt1[, list(l=unique(lower),
                      u=unique(upper),
                      minT=as.numeric(min(time)),
                      time=as.numeric(time)
    ), by=list(surv, st)]
    ### make max. time column
    dt2[, "maxT" := c(minT[2:length(minT)], NA), by=st]
    dt2 <- na.omit(dt2)
    ### add shading
    g1 <- g1 + ggplot2::geom_rect(data=dt2, ggplot2::aes(x=NULL, y=NULL,
                                       ymax=surv, ymin=l,
                                       xmax=maxT, xmin=minT,
                                       colour=st, group=st, fill=st),
                         alpha=alpha, size=fillLineSize) +
      ggplot2::geom_rect(data=dt2, ggplot2::aes(x=NULL, y=NULL,
                              ymax=u, ymin=surv,
                              xmax=maxT, xmin=minT,
                              colour=st, group=st, fill=st),
                alpha=0.05, size=fillLineSize)
  }
  ### add lines to show times where subjects censored
  if (any(dt1$n.censor >= 1)){
    g1 <- g1 + ggplot2::geom_point(data=dt1[n.censor>=1, ],
                          ggplot2::aes(x=time, y=surv),
                          shape=censShape, size=censSize)
  }
  ### palette
  ### use palette Dark2 for prominent shades
  ### (suitable for colorblind)
  ### use palette Set2 for lighter shades as large fill area
  palette <- match.arg(palette)
  if(type=="fill"){
    g1 <- g1 + ggplot2::scale_fill_brewer(type="qual", palette=palette,
                                 guide=ggplot2::guide_legend(
                                   keywidth=legSize,
                                   keyheight=legSize))
  } else {
    g1 <- g1 + ggplot2::scale_colour_brewer(type="qual", palette=palette,
                                   guide=ggplot2::guide_legend(
                                     keywidth=legSize,
                                     keyheight=legSize))
  }
  ### scales
  g1 <- g1 +
    ggplot2::scale_x_continuous(xlab, breaks = times, limits = c(0, max(object$time))) +
    ggplot2::scale_y_continuous(ylab, limits = c(0, 1)) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(title = ggplot2::element_text(size=titTextSize),
          legend.text=ggplot2::element_text(size=legTextSize),
          legend.title=ggplot2::element_text(size=legTextSize),
          axis.text = ggplot2::element_text(size = axisLabSize),
          axis.title = ggplot2::element_text(size = axisTitSize),
          legend.key = ggplot2::element_rect(colour = NA)
    )
  ### legend title
  if(type=="fill"){
    g1 <- g1 + ggplot2::labs(group=legTitle, colour=legTitle, fill=legTitle)
  } else {
    g1 <- g1 + ggplot2::labs(colour=legTitle)
  }
  ### assign color
  if(!is.null(assignCol)){
    g1 <- g1 + ggplot2::scale_color_manual(values = assignCol)
  }
  ### remove legend if required
  if(!legend) g1 <- g1 + ggplot2::theme(legend.position = c(leg.x, leg.y))
  ### p value for log-rank test (only if >=2 groups)
  if(pval & !is.null(object$strata)) {
    sd1 <- survival::survdiff(eval(object$call$formula),
                              data=eval(object$call$data))
    p1 <- stats::pchisq(sd1$chisq,
                        length(sd1$n) - 1,
                        lower.tail=FALSE)
    p1txt <- ifelse(p1 < 0.0001,
                    "Log-rank test \n p < 0.0001",
                    paste("Log-rank test \n p =", signif(p1, sigP))
    )
    g1 <- g1 + ggplot2::annotate("text",
                        x = pX * max(dt1$time),
                        y = pY,
                        label = p1txt,
                        size = ggplot2::element_text(size=legTextSize))
  }
  ### times to show on table
  tabTime <- match.arg(tabTime)
  ### use marks from existing plot
  if (tabTime=="major"){
    times1 <- ggplot2::ggplot_build(g1)$panel$ranges[[1]]$x.major_source
  } else {
    times1 <- ggplot2::ggplot_build(g1)$panel$ranges[[1]]$x.minor_source
  }
  ### data for table
  dt3 <- data.table::data.table(
    time = summary(object, times = times1, extend = TRUE)$time,
    n.risk = summary(object, times = times1, extend = TRUE)$n.risk
  )
  ### if intercept-only model
  if (is.null(object$strata)) {
    dt3[, "st" := as.factor(rep(1, length(times1)))]
  } else {
    dt3[, "st" := summary(object, times=times1, extend=TRUE)$strata]
  }
  ### change names of strata to legend labels
  if(!is.null(legLabs)) dt3[, "st" := factor(st, labels=legLabs) ]
  ### table
  ### reverse here to plot in same order as in main plot
  g2 <- ggplot2::ggplot(data=dt3, ggplot2::aes(x=time, y=rev(st), shape=rev(st))) +
    ggplot2::geom_point(size=0) +
    ggplot2::geom_text(ggplot2::aes(label=n.risk), colour=1, size=nRiskSize) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(name="Numbers at risk", limits=c(0, max(object$time)),
                       breaks=times1) +
    ### reverse here to plot in same order as in main plot
    ggplot2::scale_y_discrete(name="", breaks=as.character(levels(dt3$st)),
                     labels=rev(levels(dt3$st))) +
    #ggplot2::ggtitle(tabTitle) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10, vjust = 1), panel.grid.major = ggplot2::element_blank(),
      	panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
    		axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
    		axis.text.y = ggplot2::element_text(face = "bold", hjust = 1)) +
    #ggplot2::theme(axis.text = ggplot2::element_text(size=axisLabSize),
    #      axis.title = ggplot2::element_text(size=axisTitSize),
    #      plot.title = ggplot2::element_text(size=tabTitTextSize),
    #      legend.title = ggplot2::element_text(size=tabLegTextSize),
    #      legend.text = ggplot2::element_text(size=tabLegTextSize)
    #) +
    ggplot2::guides(shape = ggplot2::guide_legend(title=legTitle,
                                keywidht=legSize,
                                keyheight=legSize))
  ### remove legend
  if(!legend) g2 <- g2 + ggplot2::theme(legend.position = "none")
  res <- list("table"=g2,
              "plot"=g1)
  class(res) <- c("tableAndPlot", "list")
  return(res)
}


#' Put Kaplan-Meier Plot and table together using `survMisc` package in R 
#'
#' Put Kaplan-Meier Plot and table together using `survMisc` package in R.
#' 
#' @name autoplot.tableAndPlot
#' @export autoplot.tableAndPlot
#' @import data.table
#' @aliases autoplot.tableAndPlot
#' @method autoplot tableAndPlot
#' @include autoplot.R
#' @description Uses \code{gridExtra::gridArrange}
#' to arrange a plot, it's legend and a table.
#' @title Arrange and plot a survival plot, it's legend and a table.
#' 
#' @param object An object of class \code{tableAndPlot} as returned by
#' \code{autoplot.survfit}
#' @param ... Additional arguments (not implemented)
#' @param hideTabLeg Suppress table legend. If \code{supTabLeg = FALSE}
#' then a legend will also be shown for the table
#' @param plotHeight Height of plot
#' @param tabHeight Height of table
#' 
#' @return A graph, as plotted by \code{gridExtra::grid.arrange}
#' 
#' @details Arguments to \code{plotHeigth} and \code{tabHeight} are
#' best specified as fractions adding to one e.g. \eqn{0.85 + 0.15 =1}.
#' \cr \cr
#' Other \code{ggplot2::ggplot2} objects may be plotted using this
#' method. They need to be stored in a list of two (with element \code{table}
#' if \code{hideTabLeg} is to be used). The class of this list should be
#' modified with
#' \cr \cr
#' \code{class(list1) <- c("tableAndPlot", "list")}
#' @author Chris Dardis. Based on existing work by
#' R. Saccilotto, Abhijit Dasgupta, Gil Tomas and Mark Cowley.
#' @examples
#' data(colon)
#' fit <- survfit(Surv(time,status)~rx, data=colon)
#' autoplot(autoplot(fit, timeby=500, xlab="Time (Day)"))
#' mycol <- c("rx=Lev"="blue", "rx=Lev+5FU"="red", "rx=Obs"="black")
#' autoplot(autoplot(fit, timeby=500, xlab="Time (Day)", assignCol=mycol))
autoplot.tableAndPlot <- function(object, ...,
                   hideTabLeg=TRUE,
                   plotHeight=0.80,
                   tabHeight=0.20){
  stopifnot(inherits(object, "tableAndPlot"))
  if(hideTabLeg){
    object$table <- object$table +
      ggplot2::theme(legend.key.height = NULL,
             legend.key.width = NULL,
             legend.key = ggplot2::element_rect(colour = NA, fill = NA),
             legend.text = ggplot2::element_text(colour = NA),
             legend.title = ggplot2::element_text(colour = NA))
  }
  ###
  plots <- rev(object)
  grobs <- widths <- list()
  ### collect the widths for each grob of each plot
  for (i in 1:length(plots)){
    grobs[[i]] <- ggplot2::ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  ### use do.call to get the max width
  maxwidth <- do.call(grid::unit.pmax, widths)
  ### asign the max width to each grob
  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  ### plot
  do.call(gridExtra::grid.arrange, c(grobs, nrow = 2,
                                     heights=list(c(plotHeight,tabHeight))))
}