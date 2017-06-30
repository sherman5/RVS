#' probability of sharing of rare variants in a subset of families
#' @export
#'
#' @description Computing probability of sharing of rare variants in
#'  a subset of families where rare variants are seen based on precomputed
#'  family-specific rare variant sharing probabilities.
#' @details All the subsets of families of size equal or inferior to the
#'  length of not are created, and the joint probability of each such
#'  subset not sharing a rare variant and the remaining families sharing
#'  a rare variant is obtained as the product of the family-specific rare
#'  variant sharing probabilities or its complement. The function then sums
#'  the pattern probabilities inferior or equal to the probability
#'  of the observed pattern of the not families not sharing a rare variant
#'  and the remaining families sharing a rare variant.
#' @param probs sharing probabilities for all families
#' @param shared boolean vector describing if all affected subjects
#'  in the family share the variant (TRUE if all share)
#' @return P-value of the exact rare variant sharing test requiring
#'  sharing by all affected subjects
#' @examples
#'  data(samplePedigrees)
#'  notSharedFams <- c(15159, 15053, 15157)
#'  famids <- sapply(samplePedigrees, function(p) p$famid[1])
#'  shared <- famids %in% notSharedFams
#'  probs <- sapply(samplePedigrees, RVsharing)
#'  multipleFamilyPValue(probs, shared)
multipleFamilyPValue <- function(probs, shared)
{
    # check: "not" contains at least one family 
    if (sum(!shared) == 0)
        stop("number of families not sharing the RV is zero.")

    # If all families share the variant, then return 1
    if (all(shared)) return (1)
    
    # total number, and number not sharing
    nf <- length(probs)
    nnot <- sum(!shared)

    # Probability of observed data
    p.obs <- prod(probs[shared], 1 - probs[!shared])

    # Tail probability includes case where all families share the variant
    p = prod(probs)

    for (h in 1:nnot)
    {
        comb.mat <- combn(nf, h)

        # plus the cases where the probability with two families not sharing
        # the variant is less extreme than the observed
        # Compute probability for all pairs of families not sharing
        for (i in 1:ncol(comb.mat))
        {
            ptmp <- prod(probs[-comb.mat[,i]], 1 - probs[comb.mat[,i]])
            if (ptmp <= p.obs) p <- p + ptmp
        }
    }
    return(p)
}

#' depreciated function
#' @export
#' @description This function is depreciated with version >= 2.0
#'  and should not be used, instead use multipleFamilyPValue
#' @param vec a vector of names of all families where a variant is seen
#' @param not a vector of names of families where not all affected subjects
#'  share the rare variant
#' @param pshare.data a data frame with at least two of the following columns:
#'  pshare: vector of RV sharing probabilities
#'  ped.tocompute.vec: vector of names of the families whose sharing 
#'      probability is contained in pshare. The names in the arguments
#'      vec and not must be found in ped.tocompute.vec
#' @return P-value of the exact rare variant sharing test requiring
#'  sharing by all affected subjects.
get.psubset <- function(vec, not, pshare.data)
{
    warning(paste('this function is depreciated with version >= 2.0',
        'and should not be used, instead use multipleFamilyPValue'))
    names <- pshare.data$ped.tocompute.vec
    probs <- pshare.data$pshare[names %in% vec]
    probNames <- names[names %in% vec]
    shared <- !(probNames %in% not)
    return(multipleFamilyPValue(probs, shared))
}

