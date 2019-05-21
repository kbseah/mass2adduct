#' Plot mass spectrometry imaging image of a msimat object
#'
#' From an MSI dataset imported with \code{\link{msimat}}, plot an image of the
#' signal intensity for a specified peaks. The image dimensions must be
#' specified as either \code{nrow} or |code{ncol}.
#'
#' @param d msimat; MSI dataset imported with \code{\link{msimat}}
#' @param peak Name of the peak to plot. If NULL, then plot the total intensity
#'        for all peaks
#' @param nrow Number of rows in the MSI image matrix
#' @param ncol Number of columns in the MSI image matrix
#' @param col Color scale to use (passed to \code{\link{image}})
#' @param main Title of plot, otherwise the peak number
#' @param ... Other parameters to pass to \code{\link{image}}
#'
#' @export

image.msimat <- function(d,
                         peak=NULL,
                         nrow=NULL,
                         ncol=NULL,
                         col=grDevices::terrain.colors(36),
                         main=NULL,
                         ...) {
  if (is.null(nrow) & is.null(ncol)) {
    stop("No dimensions specified for the image")
  }
  if (is.null(peak)) {
    # Find the highest peak
    #idx <- which(d$peakintensities == max(d$peakintensities))
    vec <- Matrix::rowSums(d$mat)
    if (is.null(main)) { # Give image title if none already specified
      main <- "Total intensity map"
    }
  } else if (! peak %in% d$peaks) {
    stop("Specified peak not found in msimat object")
  } else {
    idx <- which(d$peaks == peak)
    vec <- d$mat[,idx] # Get data from peak vs. spot matrix
    peakname <- as.character(d$peaks[idx])
    if (is.null(main)) {
      main <- paste(c("Intensity map for peak", peakname),collapse=" ")
    }
  }
  # Check that specified dimensions fit
  if (!is.null(nrow)) {
    if (dim(d$mat)[1]%%nrow != 0) {
      stop ("Dimensions specified do not fit image, please check")
    }
    img <- matrix(vec,nrow=nrow)
  } else if (!is.null(ncol)) {
    if (dim(d$mat)[1]%%ncol != 0) {
      stop ("Dimensions specified do not fit image, please check")
    }
    img <- matrix(vec,ncol=ncol)
  }
  # Plot image
  image(img,
        col=col,
        main=main,
        ...)
}
