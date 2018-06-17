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
