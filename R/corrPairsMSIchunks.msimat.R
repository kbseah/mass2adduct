#' Calculate correlations of pairs of mass peaks from MSI data in chunks
#'
#' Calculate correlations of pairs of mass peaks from MSI data (imported to
#' R with the \code{\link{readMSI}} function). The list of pairs is supplied
#' as a data.frame with the parameter \code{pairs}, but is broken up into "chunks"
#' for processing, to avoid going over a specified memory limit.
#'
#' This function is otherwise equivalent to \code{\link{corrPairsMSI}} with
#' option how="parallel". The memory requirement is estimated conservatively but
#' has only been tested empirically (i.e. the process might go over the limit)
#' in practice.
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
#' @param ncores integer; Number of cores if using parallel version. Default is
#'        total number of cores minus one.
#' @param mem.limit integer; Memory limit estimate (in Gb). (default: 5)
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
#'         P.value - P-value for the correlation
#'
#'         Significance - Whether p-value meets cutoff (specified by "p.val", with
#'                        Bonferroni correction)
#'
#' @seealso \code{\link{corrSinglePeakMSI}} to calculate correlations of all
#'          peaks in a dataset vs. a single peak of interest
#'
#' @export

corrPairsMSIchunks.msimat <- function(d,
                                diff,
                                p.val=0.05,
                                method="pearson",
                                alternative="greater",
                                ncores=NULL,
                                mem.limit=5,
                                ...) {
    # Adapted from original code by Moritz
    # Get vectors representing parent and adduct ion masses
    ions.parent <- diff[["A"]]
    ions.adduct <- diff[["B"]]
    # Subset the MSI data frame to contain only target ions
    peaklist <- d[["peaks"]]
    idx.parent <- match(ions.parent, peaklist)
    idx.adduct <- match(ions.adduct, peaklist)
    numpairs <- length(idx.parent)

    if (is.null(ncores)) {
        # If number of cores not specified, one fewer than total detected cores
        ncores <- parallel::detectCores() - 1
    }
    
    # initialize output data frame
    df <- data.frame(Estimate=numeric(),
                     P.value=numeric())
    
    # Estimate memory required, being v conservative with 10 bytes per numeric
    # Num pixels * Num pairs * 2 * 10 B * num cores
    # However, not taking memory size of output into account...
    # Using as.numeric to avoid integer overflow
    mem.needed <- as.numeric(dim(d[["mat"]])[1]) * as.numeric(numpairs) * 2.0 * 10.0 * as.numeric(ncores) / 1e9
    numchunks <- ceiling(mem.needed/mem.limit)
    chunksize <- floor(numpairs/numchunks)
    cat (paste(c("Mem needed is ",mem.needed,"Gb","and number of chunks",numchunks,"\n"),collapse=" "))
    
    # Start chunking
    chunkidx <- 0:(numchunks-1) * chunksize
    chunkidx <- c(chunkidx, numpairs)
    for (i in 1:(length(chunkidx)-1)) {
        cat(paste(c("Processing chunk",i,"...\n"),collapse=" "))
        # Get slices of the data for this chunk
        startidx <- chunkidx[i] + 1
        stopidx <- chunkidx[i + 1]
        idx.parent.subset <- idx.parent[startidx:stopidx]
        idx.adduct.subset <- idx.adduct[startidx:stopidx]
        currchunksize <- length(idx.parent.subset)
        
        # Continue with the correlation testing
        A <- d[["mat"]][,idx.parent.subset]
        B <- d[["mat"]][,idx.adduct.subset]
        # If original matrix is a Tsparse matrix object, convert to normal
        # matrices for speed; otherwise accessing the indices at the lapply
        # function will be slow
        A <- Matrix::as.matrix(A)
        B <- Matrix::as.matrix(B)
        
        # Pairwise correlation with p-values
        cat (paste(c("Calculating correlations between",currchunksize,"pairs\n")))
    
        #cat (paste(c("Size of A is",object.size(A)/1e9,"Gb","\n"),collapse=" "))
        #cat (paste(c("Size of B is",object.size(B)/1e9,"Gb","\n"),collapse=" "))
        testresult <- parallel::mclapply(1:currchunksize,
                                         function(x) {unlist(cor.test(A[,x], B[,x], method=method, alternative=alternative, ...)[c("estimate","p.value")])},
                                         mc.cores=ncores)
        testresult <- unlist(testresult)
        testresult <- matrix(testresult,nrow=currchunksize,byrow=TRUE)
        df <- rbind(df, data.frame(testresult[,1],testresult[,2]))
        #cat (paste(c("Size of testresult is",object.size(testresult)/1e9,"Gb","\n"),collapse=" "))
    }
    
    names(df) <- c("Estimate","P.value")
    # test for significance with Bonferroni adjusted p-values (p value <= 0.05/number of pixels])
    adjusted.pvalue <- p.val/numpairs
    df$Significance <- ifelse(df$P.value<=adjusted.pvalue, 1, 0)
    # add the parental and adduct ions to the dataframe
    df$A <- ions.parent
    df$B <- ions.adduct
    df$diff <- diff[["diff"]]
    
    if (is.null(diff[["matches"]])) {
      # Reorder fields
      df <- df[,c(4,5,6,1,2,3)]
    } else {
      # If this is the output from adductMatch.massdiff, also bind the adduct names
      df$matches <- diff[["matches"]]
      # Reorder fields
      df <- df[,c(4,5,6,7,1,2,3)]
    }
    
    # Report number of significantly-correlated pairs found
    num.signif <- sum(df$Significance)
    cat ("Significant correlations found at p-value cutoff of", p.val, "(with Bonferroni correction):", num.signif, "\n")
    # Return data frame
    return(df)
}
