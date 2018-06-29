#' Summary method for class massdiff
#'
#' @param d massdiff; Object created by function \code{\link{massdiff}}
#'
#' @export

summary.massdiff <- function(d) {
  mdinfo <- data.frame(length(d[["diff"]]),
                       min(d[["A"]]),
                       max(d[["B"]]),
                       !is.null(d[["P.value"]]),
                       !is.null(d[["matches"]])
                       )
  names(mdinfo) <- c("Pairs","Lowest mass","Highest mass","Contains correlation test results","Contains matches")
  cat("Object of class massdiff:\n")
  return(mdinfo)
}
