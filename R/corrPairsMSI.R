#' Calculate correlations of pairs of mass peaks from MSI data
#'
#' Calculate correlations of pairs of mass peaks from MSI data (imported to
#' R with the \code{\link{msimat}} function). The list of pairs is supplied
#' as a data.frame with the parameter \code{pairs}.
#'
#' Example usage scenario for this function: Mass differences for all pairwise
#' combinations of masses are tabulated with \code{\link{massdiff}}, and the
#' peak pairs corresponding to a specific adduct of interest are extracted with
#' \code{\link{diffGetPeaks}}. Check which peaks are significantly correlated to
#' each other with this function.
#'
#' Correlation is calculated with \code{\link{cor.test}} function and by default
#' uses the Pearson method (two-sided, because correlations can be both positive
#' and negative). The Bonferroni correction is applied to the p-values before
#' assessing significance, but the original p-value is reported in column
#' "P.value" and can be used e.g. for false discovery rate analysis.
#'
#' @param d msimat; MSI data with peaks as columns and pixels as rows,
#'        output from \code{\link{msimat}} function
#' @param diff massdiff; List of mass differences, parent and putative adduct
#'        ions, as produced by function \code{\link{massdiff}}
#' @param p.val numeric; p-value cutoff (before Bonferroni correction) (default:
#'        0.05)
#' @param method string; Method to use for \code{\link{cor.test}} (default:
#'        "pearson")
#' @param alternative string; Which alternative for \code{\link{cor.test}}
#'        (default: "greater")
#' @param how string; How to implement the multiple correlation tests. Options:
#'        "loop" (for loop - slow), "parallel" (multiple processors, using
#'         \code{\link{mclapply}} from package \code{parallel}), or "apply"
#'        (vectorized lapply function - default). Option "parallel" uses forking
#'        processes and is therefore not suitable for Windows systems.
#' @param ncores integer; Number of cores if using \code{how="parallel"}. Default is
#'        total number of cores minus one.
#' @param ... Other parameters to pass to \code{\link{cor.test}}
#'
#' @return Object of class data.frame with the following fields:
#'
#'         A - First peak mass in each pair
#'
#'         B - Second peak mass in the pair
#'
#'         Estimate - Estimated correlation
#'
#'         P.value - P-value for the correlation (uncorrected)
#'
#'         Significance - Whether p-value meets cutoff (specified by "p.val", with
#'                        Bonferroni correction)
#'
#' @seealso \code{\link{corrPairsMSIchunks}} to split input data into manageable
#'          sizes to avoid going over available memory. Suggested if using the
#'          standard corrPairsMSI function causes R to run out of memory and
#'          crash.
#'
#' @examples
#' d <- msimat(csv=system.file("extdata","msi.csv",package="mass2adduct"),sep=";")
#' d.diff <- massdiff(d) # Calculate mass differences from imported MSI data
#' d.diff.annot <- adductMatch(d.diff,add=adducts2) # Find mass diffs matching adducts
#' d.diff.annot.cor <- corrPairsMSI(d,d.diff.annot,how="apply") # Perform correlation tests on annotated peak pairs
#'
#' @export

corrPairsMSI <- function(d,
                         diff,
                         p.val=0.05,
                         method=c("pearson","kendall","spearman"),
                         alternative=c("two.sided","less","greater"),
                         how=c("apply","parallel","loop"),
                         ncores=NULL,
                         ...) {
    if (class(d) != "msimat") {
        stop("Input parameter d must be an msimat object")
    }
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
    message("Calculating correlations between ",numpairs," pairs")

    # Run cor.test, using different strategies
    if (how[1] == "loop") { # Use a for-loop to run cor.test
        df <- data.frame(Estimate=numeric(numpairs),
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
    } else if (how[1] == "parallel") { # Use parallelized lapply to run cor.test
        if (is.null(ncores)) {
            # If number of cores not specified, one fewer than total detected cores
            ncores <- parallel::detectCores() - 1
        }
        testresult <- parallel::mclapply(1:numpairs,
                                         function(x) {unlist(cor.test(A[,x], B[,x], method=method, alternative=alternative, ...)[c("estimate","p.value")])},
                                         mc.cores=ncores)
        testresult <- unlist(testresult)
        testresult <- matrix(testresult,nrow=numpairs,byrow=TRUE)
        df <- data.frame(Estimate=testresult[,1],
                         P.value=testresult[,2])
    } else { # Default - use lapply
        testresult <- lapply(1:numpairs,
                             function(x) {unlist(cor.test(A[,x], B[,x], method=method, alternative=alternative, ...)[c("estimate","p.value")])})
        testresult <- unlist(testresult)
        testresult <- matrix(testresult,nrow=numpairs,byrow=TRUE)
        df <- data.frame(Estimate=testresult[,1],
                         P.value=testresult[,2])
    }

    # test for significance with Bonferroni adjusted p-values (p value <= 0.05/number of pixels])
    adjusted.pvalue <- p.val/numpairs
    df$Significance <- ifelse(df$P.value<=adjusted.pvalue, 1, 0)
    # add the parental and adduct ions to the dataframe
    df$A <- ions.parent
    df$B <- ions.adduct
    df$diff <- diff[["diff"]]
    if (is.null(diff[["matches"]])) {
      # Retain only relevant fields
      df <- df[,c(4,5,6,1,2,3)]
    } else {
      # If this is the output from adductMatch.massdiff, also bind the adduct names
      df$matches <- diff[["matches"]]
      # Retain only relevant fields
      df <- df[,c(4,5,6,7,1,2,3)]
    }
    # Report number of significantly correlated pairs found
    num.signif <- sum(df$Significance)
    message("Significant correlations found at p-value cutoff of ", p.val, " (with Bonferroni correction): ", num.signif)
    # Return data frame
    return(df)
}
