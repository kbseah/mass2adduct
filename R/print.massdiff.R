#' Print method for class massdiff
#'
#' @param d massdiff; Object created by function \code{\link{massdiff}}
#'
#' @export

print.massdiff <- function(d) {
    mdinfo <- data.frame(length(d[["diff"]]),
                         min(d[["A"]]),
                         max(d[["B"]]),
                         !is.null(d[["matches"]])
                         )
    names(mdinfo) <- c("Pairs","Lowest mass","Highest mass","Contains matches")
    cat("Object of class massdiff:\n")
    print(mdinfo)
    return(mdinfo)
}
