#' probability of sharing of rare variants in a subset of families
#' @export
#'
#' @description Computing probability of sharing of rare variants in
#' a subset of families where rare variants are seen based on precomputed
#' family-specific rare variant sharing probabilities.
#' @details All the subsets of families of size equal or inferior to the
#' length of not are created, and the joint probability of each such
#' subset not sharing a rare variant and the remaining families sharing
#' a rare variant is obtained as the product of the family-specific rare
#' variant sharing probabilities or its complement. The function then sums
#' the pattern probabilities inferior or equal to the probability
#' of the observed pattern of the not families not sharing a rare variant
#' and the remaining families sharing a rare variant.
#' @param sharingProbs named vector of sharing probabilties, where names
#' correspond to famid value of pedigree
#' @param observedSharing boolean vector describing if all affected subjects
#' in the family share the variant (TRUE if all share)
#' @param minPValue the minimum p-value threshold, once the true p-value is
#' determined to be less than this, the computation stops and minPValue is
#' returned - this prevents extremely long computations for extremely small
#' p-values
#' @return P-value of the exact rare variant sharing test requiring
#' sharing by all affected subjects
#' @examples
#' data(samplePedigrees)
#' probs <- sapply(samplePedigrees, RVsharing)
#' notSharedFams <- c(15159, 15053, 15157)
#' famids <- sapply(samplePedigrees, function(p) p$famid[1])
#' shared <- !famids %in% notSharedFams
#' names(shared) <- names(probs)
#' multipleFamilyPValue(probs, shared)
#' @references Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E.,
#' Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H.
#' and Ruczinski, I. (2014) Inferring rare disease risk variants based on
#' exact probabilities of sharing by multiple affected relatives.
#' Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198.
multipleFamilyPValue <- function(sharingProbs, observedSharing, minPValue=0)
{
    # check all pedigrees are given
    if (!all(names(observedSharing) %in% names(sharingProbs)))
        stop('missing some pedigrees in sharingProbs argument')

    # line up and subset vectors
    sharingProbs <- sharingProbs[names(observedSharing)]
    sharingProbs <- unname(sharingProbs)
    observedSharing <- unname(observedSharing)

    # probability of observed data
    true_pObserved <- prod(sharingProbs[observedSharing]) * 
        prod(1 - sharingProbs[!observedSharing])
    if (true_pObserved == 0)
        return(0)

    # sum probabilities of both branches of the tree
    sumBranches <- function(ndx, prod)
    {
        # get sum of probs starting at leaf of this node
        leafSum <- function(p)
        {
            # TODO account for equal sequences exactly, then check less than
            if (p < pObserved + 1e-3) return(p)
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

    # test higher observed values to see if under min p value
    pvalue <- 1
    pObserved <- 1
    while (pvalue > minPValue & pObserved > true_pObserved)
    {
        pObserved <- max(pObserved / 10, true_pObserved)
        pvalue <- sumBranches(1, 1)
    }

    if (pObserved < true_pObserved)
        stop('binary tree burn in failed')
    return(pvalue)
}

#' convert snpMatrix to a list of vectors of sharing
#' @inheritParams multipleVariantPValue
#' @return list of boolean vectors indicating sharing pattern for each variant
convertMatrix <- function(snpMat, famInfo, minorAllele=NA)
{
    mat <- snpMat@.Data
    sapply(colnames(mat), function(var)
    {
        ret <- c()
        if (is.na(minorAllele))
        {
            minorCount <- ifelse(sum(mat[,var] == 1) < sum(mat[,var] == 3), 1, 3)
        }
        else 
        {
            minorCount <- ifelse(minorAllele[var] == 1, 1, 3)
        }
    
        for (fid in unique(famInfo$pedigree))
        {
            aff_rows <- which(famInfo$pedigree == fid & famInfo$affected == 2)
            aff_alleles <- mat[aff_rows, var]
            if (!all(aff_alleles != minorCount & aff_alleles != 2))
            {
                ret <- c(ret, all(aff_alleles == minorCount | aff_alleles == 2))
                names(ret)[length(ret)] <- fid
            }
        }
        return(ret)
    }, USE.NAMES=TRUE, simplify=FALSE)
}

#' generalization of multipleFamilyPValue to multiple variants
#' @export
#'
#' @description Computes a p-value for each variant sharing pattern
#'  across families
#' @details For each variant, the families which have all sequenced subjects
#'  sharing the variant and the families which have some sequenced subjects
#'  sharing the variant are recorded. These values are passed
#'  to multipleFamilyPValue
#' @param snpMat SnpMatrix
#' @param famInfo data frame containing pedigree, member, father, mother,
#' sex, affected fields for each sequenced subject
#' @param sharingProbs vector of sharing probabilites, must be a named vector
#' with famid's for each probability
#' @param minorAllele vector specifying the minor allele of each variant
#' @param filter criteria for filtering pvalues
#' @param alpha parameter for filter
#' @return list containing p-values and potential p-values for each variant
multipleVariantPValue <- function(snpMat, famInfo, sharingProbs,
minorAllele=NA, filter=NULL, alpha=0)
{
    # convert matrix to list of families with each allele
    shareList <- convertMatrix(snpMat@.Data, famInfo, minorAllele)

    # calculate potential p-values
    pot_pvals <- sapply(shareList, function(vec)
    {
        prob <- 1
        for (fid in names(vec))
        {
            prob <- prob * unname(sharingProbs[fid])
        }
        return(prob)
    }, USE.NAMES=TRUE)

    # subset data if filter is requested
    ppval_cutoff <- 1
    if (!is.null(filter))
    {
        sorted_ppvals <- sort(unname(pot_pvals))
        cutoff <- alpha / (1:length(pot_pvals))
        ppval_cutoff <- sorted_ppvals[max(which(sorted_ppvals < cutoff))]
    }

    # compute p-values
    pvals <- sapply(names(shareList), function(var)
    {
        if (pot_pvals[var] <= ppval_cutoff)
            multipleFamilyPValue(sharingProbs, shareList[[var]])
        else
            NA
    }, USE.NAMES=TRUE)

    # return p-values and potential p-values
    return(list(pvalues=pvals[!is.na(pvals)], potential_pvalues=pot_pvals))
}

#' enrichment p-value acroos multiple families and variants
#' @export
#'
#' @description Computes a p-value for all variants seen across all families
#' @details For each variant, the families which have all sequenced subjects
#' sharing the variant and the families which have some sequenced subjects
#' sharing the variant are recorded. All unique (family, variant) pairs
#' are accumulated into a single vector and passed to multipleFamilyPValue
#' @param threshold minimum p-value threshold passed to multipleFamilyPValue
#' @inheritParams multipleVariantPValue
#' @return p-value
#' @references Fu, J., Beaty, T.H., Scott, A.F., Hetmanski, J., Parker, M.M., Bailey-Wilson, J.E.,
#' Marazita, M.L., et al. 2017. “Whole Exome Association of Rare Deletions in Multiplex Oral Cleft Families.” Genetic Epidemiology 41 (1): 61–69. doi:10.1002/gepi.22010.
enrichmentPValue <- function(snpMat, famInfo, sharingProbs, threshold=0)
{
    # convert matrix to list of families with each allele
    shareList <- convertMatrix(snpMat@.Data, famInfo)

    # make sharing events unique for each variant by renaming
    probs <- pattern <- c()
    for (var in names(shareList))
    {
        temp_probs <- sharingProbs[names(shareList[[var]])]
        temp_pattern <- shareList[[var]]
        names(temp_pattern) <- paste(var, names(temp_pattern), sep='_')
        names(temp_probs) <- paste(var, names(temp_probs), sep='_')
        probs <- c(probs, temp_probs)
        pattern <- c(pattern, temp_pattern)
    }

    # return p-value (or threshold if p-value is smaller)
    pval <- multipleFamilyPValue(probs, pattern, threshold)
    return(ifelse(pval <= threshold, threshold, pval))
}

#' deprecated function
#' @export
#' @description This function is deprecated with version >= 2.0
#' and should not be used, instead use multipleFamilyPValue
#' @param vec a vector of names of all families where a variant is seen
#' @param not a vector of names of families where not all affected subjects
#' share the rare variant
#' @param pshare.data a data frame with at least two of the following columns:
#' pshare: vector of RV sharing probabilities
#' ped.tocompute.vec: vector of names of the families whose sharing 
#' probability is contained in pshare. The names in the arguments
#' vec and not must be found in ped.tocompute.vec
#' @return P-value of the exact rare variant sharing test requiring
#' sharing by all affected subjects.
#' @examples
#' data(samplePedigrees)
#' notSharedFams <- c(15159, 15053, 15157)
#' famids <- sapply(samplePedigrees, function(p) p$famid[1])
#' notShared <- famids %in% notSharedFams
#' probs <- sapply(samplePedigrees, RVsharing)
#' get.psubset(famids, notShared, data.frame(pshare=probs,
#' ped.tocompute.vec=famids))
get.psubset <- function(vec, not, pshare.data)
{
    warning(paste('this function is deprecated',
        'and should not be used, instead use multipleFamilyPValue'))
    names.vec <- pshare.data$ped.tocompute.vec
    probs <- pshare.data$pshare[names.vec %in% vec]
    probNames <- names.vec[names.vec %in% vec]
    shared <- !(probNames %in% not)
    names(probs) <- names(shared) <- probNames
    return(multipleFamilyPValue(probs, shared))
}

