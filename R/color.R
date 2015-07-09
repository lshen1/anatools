##' diverging colour palette function with set midpoint
##'
##' returns a palette function that maps values to colours, with
##'  a midpoint (defaulting to 0) corresponding to the central colour
##' @title diverging_palette
##' @export
##' @param d data giving the range of the palette
##' @param centered logical, whether to use both sides from the midpoint symmetrically
##' @param midpoint numeric value corresponding to the central colour
##' @param colors vector of colors, length must be odd
##' @return a function
##' @examples
##' grid.raster(diverging_palette(1:10, T, mid=2, col=c("blue", "white", "red"))(1:10))
##' diverging_color(1:10, T, midpoint=2, col=c("blue", "white", "red"), breaks=100)
diverging_palette <- function(d = NULL, centered = TRUE, midpoint = 0,
                              colors = RColorBrewer::brewer.pal(7,"PRGn")){
  
  half <- length(colors)/2
  
  if(!length(colors)%%2) stop("requires odd number of colors")
  values <- if(centered) {
    low <- seq(min(d), midpoint, length=half)
    high <- seq(midpoint, max(d), length=half)
    c(low[-length(low)], midpoint, high[-1])
  } else {
    mabs <- max(abs(d - midpoint))
    seq(midpoint-mabs, midpoint + mabs, length=length(colors))
  }
  scales::gradient_n_pal(colors, values = values)
}

##' diverging colour palette function with set midpoint
##'
##' returns a palette colours that maps values to colours, with
##'  a midpoint (defaulting to 0) corresponding to the central colour
##' @title diverging_color
##' @export
##' @param d data giving the range of the palette
##' @param centered logical, whether to use both sides from the midpoint symmetrically
##' @param midpoint numeric value corresponding to the central colour
##' @param colors vector of colors, length must be odd
##' @return a set of colours that maps the values.
##' @seealso \code{\link{diverging_palette}} which this function wraps.
##' @examples
##' grid.raster(diverging_palette(1:10, T, mid=2, col=c("blue", "white", "red"))(1:10))
##' diverging_color(1:10, T, midpoint=2, col=c("blue", "white", "red"), breaks=100)
diverging_color <- function(d = NULL, centered = TRUE, midpoint = 0, breaks=10,
                              colors = RColorBrewer::brewer.pal(7,"PRGn")){
  return(c(diverging_palette(d=d, centered=centered, midpoint=midpoint,
                    colors=colors)(seq(min(d), max(d),length.out=breaks))))
}


##' Add transparent/alpha to a set of colours.
##'
##' returns a set of transparent colours based on the alpha value specified.
##' @title add_alpha
##' @export
##' @param colors vector of colors.
##' @param alpha Default is 1:No transparent. 
##' @return a set of transparent colours based on the alpha value specified.
##' @examples
##' add_alpha(colors=c("black", "red", "blue"), alpha=0.5)
add_alpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}