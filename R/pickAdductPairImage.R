#' Display parent/adduct ion pair MSI images
#'
#' Interactively display MSI images for ion pairs of interest by clicking on
#' the mass spectrum to select a specific peak (either in the role of a parent
#' or adduct ion). This allows a visual inspection of the MSI image for the
#' pair, to evaluate if the ions are spatially correlated.
#'
#' The dimensions of the plot must be known and supplied to either the
#' \code{nrow} or \code{ncol} options. Otherwise, no MSI images will be drawn
#' and only the numerical data reported. Parent and adduct ion images will be
#' plotted side by side, with a new dev object for each pair.
#'
#' @param d msimat object used to plot the mass spectrum
#' @param diff massdiff object containing the annotated mass differences of
#'        interest, with adduct matches found by \code{\link{adductMatch}}. If
#'        it has been annotated with correlation test results using
#'        \code{\link{corrPairsMSI}}, the correlation coefficient and p-value
#'        will also be used to label the plot, in addition to the adduct name
#'        and mass difference.
#' @param which character; Whether the peak being selected represents the
#'        parent or adduct ion
#' @param match character; Name of the adduct annotation to visualize. If
#'        \code{match=NULL} (default) then all annotated parent-adduct ion pairs
#'        involving this peak will be visualized (warning: there might be a lot
#'        of them!)
#' @param nrow numeric
#' @param ncol numeric; Number of rows or columns in the MS image. This will be
#'        passed to \code{\link{image.msimat}}, and at least one of them must
#'        be specified if you would like to display the images. Otherwise only
#'        the data frame will be returned, without visualization.
#'
#' @return massdiff object with all mass pairs that include the peak of interest
#'         subsetted from the object passed to argument \code{diff}
#'
#' @seealso \code{\link{pointsAdducts}} Overlay adduct or parent ion peaks
#'          from a massdiff object onto a mass spectrum plot
#' @seealso \code{\link{plot.msimat}} Plot method for msimat object to draw
#'          mass spectrum with total intensities
#' @seealso \code{\link{identify.msimat}} Identify peaks on mass spectrum by
#'          clicking on the plot
#' @seealso \code{\link{image.msimat}} Image method for msimat object, to
#'          render MSI image for peaks of interest.
#' @export

pickAdductPairImage <- function(d,
                                diff,
                                which=c("parent","adduct"),
                                match=NULL,
                                nrow=NULL,
                                ncol=NULL
                                ) {
  # Check input
  if (class(d) != "msimat") {
    stop("Argument to d should be object of class msimat")
  }
  if (!"massdiff" %in% class(diff)) {
    stop("Argument to diff should be object of class massdiff")
  }
  if (is.null(diff$match)) {
    stop("massdiff object has not been annotated with adduct matches")
  }
  if (is.null(match)) {
    message("No adduct type specified, using all")
  }

  if (!which[1] %in% c("parent","adduct")) {
    stop("Argument to which should be either \"parent\" or \"adduct\"")
  }

  # Pick the point
  peakidx <- identify.msimat(d, plot=TRUE, give.names=FALSE, n=1)

  # Subset the diff table by the desired masses
  if (which[1] == "parent") {
    diffidx <- which(diff$A == d$peaks[peakidx])
  } else if (which[1] == "adduct") {
    diffidx <- which(diff$B == d$peaks[peakidx])
  }
  diff.subset <- diff[diffidx,]
  if (is.null(match)) {
    diff.subset <- diff.subset[which(!is.na(diff.subset$matches)),]
  } else {
    if (!match %in% diff$matches) {
      stop("Adduct name not found in the massdiff object")
    }
    diff.subset <- diff.subset[which(diff.subset$matches == match),]
  }

  # Plot images if nrow or ncol supplied
  if (!is.null(nrow) | !is.null(ncol)) {
    for (i in 1:dim(diff.subset)[1]) {
      dev.new()
      par(mfrow=c(1,2))
      image(d=d,
            peak=diff.subset$A[i],
            nrow=nrow,ncol=ncol,
            main=paste(c("Parent",
                         as.character(diff.subset$A[i])),
                      collapse=" "))
      image(d=d,
            peak=diff.subset$B[i],
            nrow=nrow,ncol=ncol,
            main=paste(c("Adduct",
                         as.character(diff.subset$B[i])),
                       collapse=" "))
      matchtitle <- c("Match:", as.character(diff.subset$matches[i]), " ",
                      "Diff:", diff.subset$diff[i])
      if (!is.null(diff.subset$Estimate[i]) & !is.null(diff.subset$P.value[i])) {
        matchtitle <- c(matchtitle, " ",
                        "Estimate:", diff.subset$Estimate[i], " ",
                        "P.value:", diff.subset$P.value[i])
      }
      mtext(paste(matchtitle,
                  collapse=" "),
            side=1,
            line=-2,
            outer=TRUE,
            cex=0.75)
    }
  }

  # Return massdiff subset
  return(diff.subset)
}
