#' Correlation-weighted histogram of massdiff object
#'
#' Bin mass differences from a massdiff object produced by
#' \code{\link{massdiff}} into a histogram, with counts weighted by the
#' correlation estimate. The massdiff object must contain the correlation test
#' results, product of \code{\link{corrPairsMSI}} or
#' \code{\link{corrPairsMSIchunks}} functions.
#'
#' @param diff massdiff; Output from function \code{\link{massdiff}},
#'        containing three numeric vectors (A, B, diff) representing mass
#'        peaks and their differences, processed with \code{\link{corrPairsMSI}}
#'        which adds columns Estimate, P.value, and Significance
#' @param widthFunc character; function to use to bin mass differences into
#'        histogram. (default "equal", other options to be added)
#' @param width numeric; bin width in m/z units (default 0.01)
#' @param ncores integer; number of processors to use for parallelized
#'        implementation. If NULL (default) then run with single processor
#' @param ... Other options to be passed to hist()
#'
#' @return Object of classes hist and massdiffhist
#' @seealso \code{\link{massdiff}} to generate the mass difference list
#' @seealso \code{\link{plot.massdiffhist}} for plotting and annotating the
#'          resulting object
#' @export

histWeighted <- function(diff, widthFunc="equal", width=0.01, ncores=NULL, ...) {
  if (! "massdiff" %in% class(diff)) {
    stop("Input must be of class massdiff")
  }
  if (is.null(diff$Estimate)) {
    stop("Requires massdiff object with correlation estimates. Use hist() for plain massdiff")
  }
  if (!is.numeric(diff$diff)) {
    stop("Input mass difference list must be numeric")
  }
  # Make a standard massdiff object
  output <- hist(diff, widthFunc=widthFunc, width=width, ...)
  weightedcounts <- vector(mode="numeric",length=length(output$counts))
  
  if (!is.null(ncores)) {
    # Parallelized implementation
    #ncores <- parallel::detectCores() - 1
    weightedcounts <- parallel::mclapply(1:length(weightedcounts),
                                         function(x) { sum(diff$Estimate[which(diff$diff > output$breaks[x]
                                                                               & diff$diff <= output$breaks[x+1])])},
                                         mc.cores=ncores)
  } else {
    weightedcounts <- lapply(1:length(weightedcounts),
                             function(x) { sum(diff$Estimate[which(diff$diff > output$breaks[x]
                                                                   & diff$diff <= output$breaks[x+1])])})
  }
  
  weightedcounts <- unlist(weightedcounts)
  weighteddensity <- weightedcounts / sum(weightedcounts)
  
  #return(weightedcounts)
  
  output$counts <- weightedcounts
  output$density <- weighteddensity
  
  class(output) <- c("massdiffhist","histogram")
  return(output)
}
