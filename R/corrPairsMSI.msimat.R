#' Calculate correlations of pairs of mass peaks from MSI data
#'
#' Calculate correlations of pairs of mass peaks from MSI data (imported to
#' R with the \code{\link{readMSI}} function). The list of pairs is supplied
#' as a data.frame with the parameter \code{pairs}.
#'
#' Example usage scenario for this function: Mass differences for all pairwise
#' combinations of masses are tabulated with \code{\link{diffTabulate}}, and the
#' peak pairs corresponding to a specific adduct of interest are extracted with
#' \code{\link{diffGetPeaks}}. Check which peaks are significantly correlated to
#' each other with this function.
#'
#' Correlation is calculated with \code{\link{cor.test}} function from the
#' \code{stats} package, and by default uses the Pearson methone (one-sided, i.e.
#' positive correlations). The Bonferroni correction is applied to the p-values
#' before assessing significance, but the original p-value is reported in column
#' P.value.
#'
#' @param d msimat; MSI data with peaks as columns and pixels as rows,
#'        output from \code{\link{readMSI}} function
#' @param diff massdiff; List of mass differences, parent and putative adduct
#'        ions, as produced by function \code{\link{diffTabulate}}
#' @param p.val numeric; p-value cutoff (before Bonferroni correction) (default: 0.05)
#' @param method string; Method to use for \code{\link{cor.test}} (default: "pearson")
#' @param alternative string; Which alternative for \code{\link{cor.test}} (default: "greater")
#' @param how string; How to implement the multiple correlation tests. Options:
#'        "loop" (for loop - slow), "parallel" (multiple processors, using
#'         \code{\link{mclapply}} from package \code{parallel}), or "apply"
#'        (vectorized lapply function - default).
#' @param ncores integer; Number of cores if using parallel version. Default is
#'        total number of cores minus one.
#' @param ... Other parameters to pass to \code{\link{cor.test}}
#'
#' @return Object of class data.frame with the following fields:
#'
#'         label - Label for row
#'
#'         A - First peak mass in each pair
#'
#'         B - Second peak mass in the pair
#'
#'         Estimate - Estimated correlation
#'
#'         P.value - P-value for the correlation
#'
#'         Significance - Whether p-value meets cutoff (specified by "p.val", with
#'                        Bonferroni correction)
#'
#' @seealso \code{\link{corrSinglePeakMSI}} to calculate correlations of all
#'          peaks in a dataset vs. a single peak of interest
#'
#' @export

corrPairsMSI.msimat <- function(d,
                                diff,
                                p.val=0.05,
                                method="pearson",
                                alternative="greater",
                                how="apply",
                                ncores=NULL,
                                ...) {
    # Adapted from original code by Moritz
    # Get vectors representing parent and adduct ion masses
    ions.parent <- diff[["A"]]
    ions.adduct <- diff[["B"]]
    # Subset the MSI data frame to contain only target ions
    peaklist <- d[["peaks"]]
    idx.parent <- match(ions.parent, peaklist)
    idx.adduct <- match(ions.adduct, peaklist)
    A <- d[["mat"]][,idx.parent]
    B <- d[["mat"]][,idx.adduct]
    # Convert DFs to normal matrices for speed; otherwise if it is a sparse
    # matrix object, accessing the indices at the lapply function will be slow
    A <- Matrix::as.matrix(A)
    B <- Matrix::as.matrix(B)
    # Pairwise correlation with p-values
    numpairs <- dim(B)[2]
    cat (paste(c("Calculating correlations between",numpairs,"pairs\n")))

    # Run cor.test, using different strategies
    if (how == "loop") { # Use a for-loop to run cor.test
        df <- data.frame(label=paste("Ion pair", 1:numpairs),
                         Estimate=numeric(numpairs),
                         P.value=numeric(numpairs),
                         Significance=numeric(numpairs))
        for (i in 1:numpairs){
            test <- cor.test(A[,i],
                             B[,i],
                             method=method,
                             alternative=alternative,
                             ...)
            df$Estimate[i] <- test$estimate
            df$P.value[i] <- test$p.value
        }
    } else if (how == "parallel") { # Use parallelized lapply to run cor.test
        if (is.null(ncores)) {
            # If number of cores not specified, one fewer than total detected cores
            ncores <- parallel::detectCores() - 1
        }
        df <- data.frame(label=paste("Ion pair", 1:numpairs))
        testresult <- parallel::mclapply(1:numpairs,
                                         function(x) {test <- cor.test(A[,x], B[,x], method=method, alternative=alternative, ...)
                                                      return(c(test$estimate,test$p.value)) },
                                         mc.cores=ncores)
        testresult <- unlist(testresult)
        testresult <- matrix(testresult,nrow=numpairs,byrow=TRUE)
        df$Estimate <- testresult[,1]
        df$P.value <- testresult[,2]
        rm(testresult)
    } else { # Default - use lapply
        df <- data.frame(label=paste("Ion pair", 1:numpairs))
        testresult <- lapply(1:numpairs,
                             function(x) {test <- cor.test(A[,x],
                                                           B[,x],
                                                           method=method,
                                                           alternative=alternative,
                                                           ...)
                                          return(c(test$estimate,test$p.value))
                                          })
        testresult <- unlist(testresult)
        testresult <- matrix(testresult,nrow=numpairs,byrow=TRUE)
        df$Estimate <- testresult[,1]
        df$P.value <- testresult[,2]
        rm(testresult)
    }

    # test for significance with Bonferroni adjusted p-values (p value <= 0.05/number of pixels])
    adjusted.pvalue <- p.val/numpairs
    df$Significance <- ifelse(df$P.value<=adjusted.pvalue, 1, 0)
    # add the parental and adduct ions to the dataframe
    df$A <- ions.parent
    df$B <- ions.adduct
    # Retain only relevant fields
    df <- df[,c(1,5,6,2,3,4)]
    # Report number of significantly correlated pairs found
    num.signif <- sum(df$Significance)
    cat ("Significant correlations found at p-value cutoff of", p.val, "(with Bonferroni correction):", num.signif, "\n")
    # Return data frame
    return(df)
}
