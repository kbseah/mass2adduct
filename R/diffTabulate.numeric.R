#' @export

diffTabulate.numeric <- function(d) {
    # Check that d is numeric
    if (!is.numeric(d)) {
        cat ("Error: Input should be numeric\n")
    } else {
        len <- length(d)
        cat ("Input is numeric vector of length", len, "with", choose(len,2), "possible pairs\n")
        if (len > 65536) { # Warning if vector is too long
            cat ("Warning: Input is longer than 65536 elements\n")
            cat ("... perhaps you should filter your peak list\n")
        } else {
            y <- comb2M(d) # Alternative to R built-in combn function
            # Format results as data frame with column labels
            y <- as.data.frame (t(y))
            names(y) <- c("A","B","diff")
            class(y) <- "massdiff"
            return(y)
        }
    }
}
