#' Print method for class msimat
#'
#' @param d msimat; Object created by function \code{\link{readMSI}}
#'
#' @export

print.msimat <- function(d) {
    matinfo <- data.frame(length(d[["peaks"]]),
                          min(d[["peaks"]]),
                          max(d[["peaks"]]),
                          length(d[["spots"]])
                          )
    names(matinfo) <- c("Peaks","Lowest mass","Highest mass","Pixels")
    cat("Object of class msimat:\n")
    print(matinfo)
    return(matinfo)
}
