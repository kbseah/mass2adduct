#' Tabulate pairwise combinations of mass peaks from imported MSI data
#'
#' Takes all pairwise combinations of mass peaks from MS imaging data imported
#' with the function \code{\link{readMSI}}.
#'
#' The number of possible pairs grows combinatorially and the calculation will
#' fail for >65536 values (the function uses a 16-bit integer internally).
#'
#' @param d numeric vector representing m/z values of mass peaks
#'
#' @return data.frame of all pairs of masses and their respective differences
#' @seealso \code{\link{diffHist}} to plot histogram of mass differences
#' @export

diffTabulate.msimat <- function(d) {
    # Extract peaklist from msimat object and apply diffTabulate.numeric
    peaklist <- d[["peaks"]]
    if (is.numeric(peaklist)) {
        output <- diffTabulate.numeric(d=peaklist)
        return(output)
    }
    else {
        cat ("Error: Input does not seem to be correct\n")
    }
}
