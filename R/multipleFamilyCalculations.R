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
    if (length(names(observedSharing)) == 0 | length(names(sharingProbs)) == 0)
        stop('sharingProbs and observedSharing must be named vectors')
    if (!all(names(observedSharing) %in% names(sharingProbs)))
        stop('missing some pedigrees in sharingProbs argument')
    # line up and subset vectors
    sharingProbs <- sharingProbs[names(observedSharing)]
    sharingProbs <- unname(sharingProbs)
    observedSharing <- unname(observedSharing)
    # calculate p-value with appropiate backend
    multipleFamilyPValue_R_Backend(sharingProbs, observedSharing, minPValue)
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
minorAllele=NULL, filter=NULL, alpha=0)
{
    # subset snpMat to only affected subjects
    if (nrow(snpMat) != length(famInfo$affected))
        stop("dimension mismatch: snpMat and famInfo")
    #snpMat <- snpMat[famInfo$affected == 2,]
    #famInfo <- famInfo[famInfo$affected == 2,]
    # check inputs
    if (!is.null(minorAllele) & length(minorAllele) != ncol(snpMat))
    {
        stop("Mismatch between number of minor alleles and size of the SnpMatrix")
    }
    multipleVariantPValue_R_Backend(snpMat, as.character(famInfo$pedigree), sharingProbs,
        minorAllele, filter, alpha)
}

#' enrichment p-value across multiple families and variants
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
#' @references Fu, J., Beaty, T.H., Scott, A.F., Hetmanski, J., Parker, M.M.,
#' Bailey-Wilson, J.E., Marazita, M.L., et al. 2017. Whole Exome Association of
#' Rare Deletions in Multiplex Oral Cleft Families. Genetic Epidemiology 41
#' (1): 61–69. doi:10.1002/gepi.22010.
enrichmentPValue <- function(snpMat, famInfo, sharingProbs, threshold=0)
{
    # determine the minor allele for each variant
    minorAllele <- sapply(colnames(snpMat), function(var)
    {
        ifelse(sum(snpMat[,var] == 1) < sum(snpMat[,var] == 3), 1, 3)
    })
    validVariants <- sapply(colnames(snpMat), function(var)
    {
        sum(snpMat[,var] == minorAllele[var] | snpMat[,var] == 2) > 0
    })
    message(paste0("Ignoring ", sum(!validVariants), " variants not present in any subject"))
    snpMat <- snpMat[,validVariants]
    minorAllele <- minorAllele[validVariants]

    # subset snpMat to only affected subjects
    if (nrow(snpMat) != length(famInfo$affected))
        stop("dimension mismatch: snpMat and famInfo")
    snpMat <- snpMat[famInfo$affected == 2,]

    enrichmentPValue_R_Backend(snpMat@.Data, as.character(famInfo$pedigree), sharingProbs,
        minorAllele, threshold)
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
