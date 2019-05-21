#' Plot massdiffhist histogram with named adducts labeled
#'
#' Plot massdiff histogram with option of labeling peaks corresponding to known
#' adducts or chemical transformations. The best-matching adducts for each
#' histogram bin are found internally with \code{\link{adductMatch.histogram}}.
#'
#' @param hist Histogram of a massdiff object produced by
#'        \code{\link{hist.massdiff}}
#' @param add data.frame; adduct data to use to annotate the mass diff peaks,
#'        (default: built-in data set \code{adducts})
#' @param labels numeric; label histogram peaks with top N adducts from the
#'        adduct table specified to option \code{add} (default: NULL - no labels)
#' @param pch Plot character to use for highlighting adduct peaks
#' @param col Color for labels and text
#' @param cex Character scaling for labelled points and text
#' @param pos Position of text label, values 1, 2, 3, 4 are below, left, above,
#'        and to the right of points respectively.
#' @param main character; Title for histogram plot
#' @param xlab character; X-axis label
#' @param ... Other parameters to pass to plot
#'
#' @export

plot.massdiffhist <- function(hist,
                              add=mass2adduct::adducts,
                              labels=NULL,
                              pch=0,
                              col="blue",
                              cex=0.5,
                              pos=4,
                              main="Mass difference histogram",
                              xlab="mass difference (m/z)",
                              ...) {
  # Plot histogram
  graphics:::plot.histogram(hist, main=main, xlab=xlab, ...)
  # Overlay points and text labels
  # Calculate top adducts
  if (!is.null(labels)) {
    matches <- adductMatch(x=hist, add=add, density=FALSE)
    matches <- matches[order(matches$counts,decreasing=T),]
    matches <- matches[1:labels,] # Take the top N hits
    points(x=matches$mass,
          y=matches$counts,
          pch=pch,
          col=col,
          cex=cex)
    graphics::text(x=matches$mass,
                   y=matches$counts,
                   label=as.character(matches$name),
                   pos=pos,
                   col=col,
                   cex=cex)
  }
}
