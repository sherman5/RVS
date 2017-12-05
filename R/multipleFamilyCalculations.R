#' probability of sharing of rare variants in a subset of families
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
#' @param sharingProbs sharing probabilities for all families
#' @param observedSharing boolean vector describing if all affected subjects
#' in the family share the variant (TRUE if all share)
#' @return P-value of the exact rare variant sharing test requiring
#' sharing by all affected subjects
#' @examples
#' data(samplePedigrees)
#' notSharedFams <- c(15159, 15053, 15157)
#' famids <- sapply(samplePedigrees, function(p) p$famid[1])
#' shared <- !famids %in% notSharedFams
#' probs <- sapply(samplePedigrees, RVsharing)
#' multipleFamilyPValue(probs, shared)
#' @references Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E.,
#' Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H.
#' and Ruczinski, I. (2014) Inferring rare disease risk variants based on
#' exact probabilities of sharing by multiple affected relatives.
#' Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198.
binaryTreePValue <- function(sharingProbs, observedSharing, minPValue=0)
{
    # check vector lengths
    if (length(sharingProbs) != length(observedSharing))
        stop('sharing pattern different length than sharing probs')

    # remove name from sharing prob
    sharingProbs <- unname(sharingProbs)

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

#' Computes a p-value for a single variant across multiple families
#' @export
#'
#' @description Computes a p-value for a single variant sharing pattern
#'  across families
#' @param sharingProbs vector of sharing probabilties
#' @param observedSharing vector of T/F indicating which pedigrees have the
#'  variant seen in all carriers
#' @param pedList list of pedigree objects
#' @return vector : p-value for each variant sharing pattern
multipleFamilyPValue <- function(sharingProbs=NULL, observedSharing, pedList=NULL)
{
    if (is.null(sharingProbs) & is.null(pedList))
        stop('must provide either sharing probilities or list of pedigrees')

    if (!is.null(sharingProbs))
        binaryTreePValue(sharingProbs, observedSharing)
    else
        binaryTreePValue(suppressMessages(RVsharing(pedList)), observedSharing)
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
#' @param share_matrix matrix entry [i,j] is 1 if subject i has
#'  variant j, otherwise it is zero
#' @param famid vector of family id's - same length as rows in share_matrix
#' @param subjects vector of subject id's - same length as rows in share_matrix
#' @param peds list of pedigrees referenced in famid
#' @return vector : p-value for each variant sharing pattern
multipleVariantPValue <- function(share_matrix, famid, subjects, peds)
{
    # check arguments
    if (length(unique(famids)) != length(peds))
        stop('need same number of pedigrees and unique family ids')

    # initialize return vector
    pvals <- rep(NA, ncol(share_matrix))

    # get family info
    unique_famid <- unique(famid)
    num_subjects <- sapply(unique_famid,  function(id)
        length(unique(subjects[famid==id])))

    # locate rows at which a new family starts
    family_start_rows <- c()
    r <- 1
    while (r <= length(famid))
    {
        family_start_rows <- c(family_start_rows, r)
        r <- r + num_subjects[which(unique_famid == famid[r])]
    }

    # calculate p-value for each variant
    for (var in 1:ncol(share_matrix))
    {
        share_probs <- c()
        full_share <- c()
        for (fam in family_start_rows)    
        {
            # get info about the carriers and variants
            sz <- num_subjects[which(unique_famid == famid[fam])]
            variants <- share_matrix[fam:(fam+sz-1), var]
            fam_num <- which(sapply(peds, function(p) p$famid[1]==famid[fam]))
            carriers <- subjects[fam:(fam+sz-1)]
            carrier_ids <- which(peds[[fam_num]]$id %in% carriers)

            # set the affected attribute for this pedigree
            peds[[fam_num]]$affected <- rep(0, length(peds[[fam_num]]$affected))
            peds[[fam_num]]$affected[carrier_ids] <- 1

            # check if all subjects share the variant or not
            if (sum(variants) == sz)
            {
                share_probs <- c(share_probs, suppressMessages(RVsharing(
                    peds[[fam_num]], useAffected=TRUE)))
                full_share <- c(full_share, TRUE)
            }
            else if (sum(variants) != 0)
            {
                share_probs <- c(share_probs, suppressMessages(RVsharing(
                    peds[[fam_num]], useAffected=TRUE)))
                full_share <- c(full_share, FALSE)
            }
        }

        # make sure at least one family has the variant
        if (length(share_probs) == 0)
            pvals[var] <- 1
        else
            pvals[var] <- binaryTreePValue(share_probs, full_share)
    }

    # return vector of p-values
    names(pvals) <- colnames(share_matrix)
    return(pvals)
}

#' enrichment p-value acroos multiple families and variants
#' @export
#'
#' @description Computes a p-value for all variants seen across all families
#' @details For each variant, the families which have all sequenced subjects
#'  sharing the variant and the families which have some sequenced subjects
#'  sharing the variant are recorded. All unique (family, variant) pairs
#'  are accumulated into a single vector and passed to multipleFamilyPValue
#' @param share_matrix matrix entry [i,j] is 1 if subject i has
#'  variant j, otherwise it is zero
#' @param famid vector of family id's - same length as rows in share_matrix
#' @param subjects vector of subject id's - same length as rows in share_matrix
#' @param peds list of pedigrees referenced in famid
#' @return vector : p-value for each variant sharing pattern
enrichmentPValue <- function(share_matrix, famid, subjects, peds, threshold)
{
    # get family info
    unique_famid <- unique(famid)
    num_subjects <- sapply(unique_famid,  function(id)
        length(unique(subjects[famid==id])))

    # locate rows at which a new family starts
    family_start_rows <- c()
    r <- 1
    while (r <= length(famid))
    {
        family_start_rows <- c(family_start_rows, r)
        r <- r + num_subjects[which(unique_famid == famid[r])]
    }

    # vectors for all (family, variant) pairs seen
    all_share_probs <- c()
    all_full_share <- c()

    # calculate p-value for each variant
    for (var in 1:ncol(share_matrix))
    {
        share_probs <- c()
        full_share <- c()
        for (fam in family_start_rows)    
        {
            # get info about the carriers and variants
            sz <- num_subjects[which(unique_famid == famid[fam])]
            variants <- share_matrix[fam:(fam+sz-1), var]
            fam_num <- which(sapply(peds, function(p) p$famid[1]==famid[fam]))
            carriers <- subjects[fam:(fam+sz-1)]
            carrier_ids <- which(peds[[fam_num]]$id %in% carriers)

            # set the affected attribute for this pedigree
            peds[[fam_num]]$affected <- rep(0, length(peds[[fam_num]]$affected))
            peds[[fam_num]]$affected[carrier_ids] <- 1

            # check if all subjects share the variant or not
            if (sum(variants) == sz)
            {
                share_probs <- c(share_probs, suppressMessages(RVsharing(
                    peds[[fam_num]], useAffected=TRUE)))
                full_share <- c(full_share, TRUE)
            }
            else if (sum(variants) != 0)
            {
                share_probs <- c(share_probs, suppressMessages(RVsharing(
                    peds[[fam_num]], useAffected=TRUE)))
                full_share <- c(full_share, FALSE)
            }
        }

        all_share_probs <- c(all_share_probs, share_probs)
        all_full_share <- c(all_full_share, full_share)
    }

    # return vector of p-values
    pval <- binaryTreePValue(all_share_probs, all_full_share, threshold)
    if (pval <= threshold)
        pval <- threshold
    return(pval)
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
    names <- pshare.data$ped.tocompute.vec
    probs <- pshare.data$pshare[names %in% vec]
    probNames <- names[names %in% vec]
    shared <- !(probNames %in% not)
    return(binaryTreePValue(probs, shared))
}

