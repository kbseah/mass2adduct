#' Calculate correlations of single target ion to other mass peaks from MSI data
#'
#' Calculate correlations of a target ion mass to all other mass peaks from
#' MSI data (imported to R with the \code{\link{readMSI}} function). The target
#' mass is supplied with the parameter \code{peak} and must match one of the
#' mass values in the MSI dataset supplied with parameter \code{d}.
#'
#' Correlation is calculated with \code{\link{cor.test}} function from the
#' \code{stats} package, and uses the default two-sided Pearson method. The
#' Bonferroni correction is applied to the p-values before assessing significance,
#' but the original p-value is reported in column P.value.
#'
#' @param d data.frame; MSI data with peaks as columns and pixels as rows,
#'        output from readMSI() function
#' @param peak numeric; Target ion for which to calculate correlations vs other ions
#'        in MSI data (input d). 
#' @param p.val numeric; p-value cutoff (before Bonferroni correction) (default: 0.05)
#' @param ... Other parameters to pass to \code{cor.test}
#'
#' @return Object of class data.frame with the following fields:
#'
#'         label - Label for row
#' 
#'         ion_of_interest - Target ion peak mass specified with parameter "peak"
#'
#'         all_ions - Peak masses for other ions in MSI data "d"
#'
#'         Estimate - Estimated correlation
#' 
#'         P.value - P-value for the correlation (without Bonferroni correction)
#'
#'         Significance - Whether p-value meets cutoff (specified by "p.val", with
#'                        Bonferroni correction)
#'
#' @seealso \code{\link{corrPairsMSI}} to calculate correlations of pairs of ions
#'          e.g. the output from \code{\link{diffTabulate}}
#'
#' @export

corrSinglePeakMSI <- function(d, peak=NULL, p.val=0.05, ...) {
    # Adapted from code by Moritz
    if (is.null(peak)) {
        # If no value specified
        cat ("Error: Ion peak mass not specified\n")
    } else {
        if (! peak %in% names(d)) {
            # If peak not in the MSI data set
            cat ("Error: Ion peak mass does not match any peak name in MSI data\n")
        } else {
            # Subset data to only ion of interest
            idx <- match(peak, names(d))
            d.peak <- d[,idx]
            # Names of all ions
            all.ions <- names(d)
            # Convert data frames to matrices
            A <- as.matrix (d.peak)
            B <- as.matrix (d)
            # Pairwise correlation
            #pairwise.corr <- sapply(seq.int(dim(A)[2]),
            #                        function(i) cor(A[,i], B[,i]))
            numpairs <- dim(B)[2]
            df <- data.frame(label=paste("Ion pair", 1:numpairs),
                             Estimate=numeric(numpairs),
                             P.value=numeric(numpairs),
                             Significance=numeric(numpairs))
            for (i in 1:numpairs){
                test <- cor.test(A[,1], B[,i], ...)
                df$Estimate[i] = test$estimate
                df$P.value[i] = test$p.value
            }
            # test for significance (p value <= 0.05/number of pixels]) 1 meaning significance, 0 meaning no significance
            adjusted.pvalue <- p.val/numpairs
            df$Significance <- ifelse(df$P.value<=adjusted.pvalue, 1, 0)
            # add the parental and adduct ions to the dataframe
            df["Ion of interest"] = peak
            df["all ions"] = all.ions
            # Extract only relevant columns
            df=df[,c(1,5,6,2,3,4)]
            # Report number of significantly correlated pairs found
            num.signif <- sum(df$Significance)
            cat ("Significant correlations found at p-value cutoff of", p.val, "(with Bonferroni correction):", num.signif, "\n")
            # Return data frame
            return(df)
        }
    }
}


