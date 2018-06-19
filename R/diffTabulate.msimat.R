#' @export

diffTabulate.msimat <- function(d) {
    # Extract peaklist from msimat object and apply diffTabulate.numeric
    peaklist <- d[["peaks"]]
    if (!is.numeric(peaklist)) {
        stop("Input peak list is not numeric\n")
    }
    output <- diffTabulate.numeric(d=peaklist)
    return(output)
}
