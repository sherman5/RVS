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
#' @param sharingProbs sharing probabilities for all families
#' @param observedSharing boolean vector describing if all affected subjects
#'  in the family share the variant (TRUE if all share)
#' @return P-value of the exact rare variant sharing test requiring
#'  sharing by all affected subjects
#' @examples
#'  data(samplePedigrees)
#'  notSharedFams <- c(15159, 15053, 15157)
#'  famids <- sapply(samplePedigrees, function(p) p$famid[1])
#'  shared <- !famids %in% notSharedFams
#'  probs <- sapply(samplePedigrees, RVsharing)
#'  multipleFamilyPValue(probs, shared)
#' @references Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E.,
#'  Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H.
#'  and Ruczinski, I. (2014) Inferring rare disease risk variants based on
#'  exact probabilities of sharing by multiple affected relatives.
#'  Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198.
multipleFamilyPValue <- function(sharingProbs, observedSharing)
{
    # remove name from sharing prob
    sharingProbs <- unname(sharingProbs)

    # probability of observed data
    pObserved <- prod(sharingProbs[observedSharing]) * 
        prod(1 - sharingProbs[!observedSharing])

    # sum probabilities of both branches of the tree
    sumBranches <- function(ndx, prod)
    {
        # get sum of probs starting at leaf of this node
        leafSum <- function(p)
        {
            if (p <= pObserved) return(p)
            else                return(sumBranches(ndx + 1, p))
        }
    
        # base case at end of tree == this path not extreme
        if (ndx > length(sharingProbs)) return(0)
        
        # multiply cumulative prob by each marginal prob
        prodLeft <- sharingProbs[ndx] * prod
        prodRight <- (1 - sharingProbs[ndx]) * prod

        # return sum of both directions
        return(leafSum(prodLeft) + leafSum(prodRight))
    }
    
    # start at root with probability 1
    return(sumBranches(1, 1))
}

#' deprecated function
#' @export
#' @description This function is deprecated with version >= 2.0
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
#' @examples
#'  data(samplePedigrees)
#'  notSharedFams <- c(15159, 15053, 15157)
#'  famids <- sapply(samplePedigrees, function(p) p$famid[1])
#'  notShared <- famids %in% notSharedFams
#'  probs <- sapply(samplePedigrees, RVsharing)
#'  get.psubset(famids, notShared, data.frame(pshare=probs, ped.tocompute.vec=famids))
get.psubset <- function(vec, not, pshare.data)
{
    warning(paste('this function is deprecated with version >= 2.0',
        'and should not be used, instead use multipleFamilyPValue'))
    names <- pshare.data$ped.tocompute.vec
    probs <- pshare.data$pshare[names %in% vec]
    probNames <- names[names %in% vec]
    shared <- !(probNames %in% not)
    return(multipleFamilyPValue(probs, shared))
}

