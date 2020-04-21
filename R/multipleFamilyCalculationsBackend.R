#### Helper functions for using a list as a stack

stackCreate <- function()
{
    stack <- list()
    stack[[1]] <- 0
    return(stack)
}

stackPushBack <- function(stack, element)
{
    size <- stack[[1]]
    stack[[size + 2]] <- element
    stack[[1]] <- size + 1
    return(stack)
}

stackGetBack <- function(stack)
{
    size <- stack[[1]]
    return(stack[[size + 1]])
}

stackPopBack <- function(stack)
{
    size <- stack[[1]]
    if (size == 0)
    {
        stop("Tried removing from empty stack")
    }
    stack[[1]] <- size - 1
    return(stack)
}

stackIsEmpty <- function(stack)
{
    size <- stack[[1]]
    return(size == 0)
}

#### Backend Functions

sumBranches <- function(index, product, allSharingProbs, observedSharingProb)
{
    pvalue <- 0
    stack <- stackCreate()
    stack <- stackPushBack(stack, c(index, product))
    while (!stackIsEmpty(stack))
    {
        current <- stackGetBack(stack)
        stack <- stackPopBack(stack)
        if (current[1] <= length(allSharingProbs))
        {
            leftProduct <- allSharingProbs[current[1]] * current[2]
            rightProduct <- (1 - allSharingProbs[current[1]]) * current[2]
            for (product in c(leftProduct, rightProduct))
            {
                if (product < observedSharingProb + .Machine$double.eps)
                {
                    pvalue <- pvalue + product
                }
                else
                {
                    stack <- stackPushBack(stack, c(current[1] + 1, product))
                }
            }
        }
    }
    return(pvalue)
}

#' R backend for multipleFamilyPValue calculation
#' @inheritParams multipleFamilyPValue
#' @return p-value
multipleFamilyPValue_R_Backend<- function(sharingProbs, observedSharing, minPValue=0)
{
    # probability of observed data
    trueObservedSharingProb <- prod(sharingProbs[observedSharing]) * 
        prod(1 - sharingProbs[!observedSharing])
    if (trueObservedSharingProb == 0)
    {
        warning("Observed sharing sequence has zero probability")
        return(0)
    }
    # calculate the exact p-value by default
    if (minPValue == 0)
    {
        return(sumBranches(1, 1, sharingProbs, trueObservedSharingProb))
    }
    # if a minimum p-value threshold is specified, start at a high observed sharing
    # probability and descend until either the p-value thresold is crossed or 
    # the true observed sharing probability is reached. This prevents spending large
    # amounts of compute time for extremely small p-values.
    pvalue <- 1
    observedSharingProb <- 1
    while (pvalue > minPValue & observedSharingProb > trueObservedSharingProb)
    {
        observedSharingProb <- max(observedSharingProb / 10, trueObservedSharingProb)
        pvalue <- sumBranches(1, 1, sharingProbs, observedSharingProb)
    }
    return(pvalue)
}

#' R backend for multipleVariantPValue calculation
#' @inheritParams multipleVariantPValue
#' @param famIds family ids corresponding to rows of the snpMap
#' @return list of p-values and potential p-values
multipleVariantPValue_R_Backend <- function(snpMat, famIds, sharingProbs,
minorAllele, filter=NULL, alpha=0)
{
    # check inputs
    if (is.null(names(sharingProbs)))
    {
        stop('sharingProbs must be a named vector')
    }
    if (is.null(minorAllele))
    {
        minorAllele <- sapply(colnames(snpMat), function(var)
        {
            ifelse(sum(snpMat[,var] == 1) < sum(snpMat[,var] == 3), 1, 3)
        })
    }
    validVariants <- sapply(colnames(snpMat), function(var)
    {
        sum(snpMat[,var] == minorAllele[var] | snpMat[,var] == 2) > 0
    })
    message(paste0("Ignoring ", sum(!validVariants), " variants not present in any subject"))
    snpMat <- snpMat[,validVariants]
    minorAllele <- minorAllele[validVariants]
    # convert matrix to list of families with each allele
    shareList <- convertMatrix(snpMat@.Data, famIds, minorAllele)
    # calculate potential p-values
    pot_pvals <- sapply(shareList, function(vec)
    {
        if (!all(names(vec) %in% names(sharingProbs)))
        {
            stop('sharingProbs is missing a value for some families')
        }
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
        {
            multipleFamilyPValue(sharingProbs, shareList[[var]])
        }
        else
        {
            NA
        }
    }, USE.NAMES=TRUE)
    # return p-values and potential p-values
    return(list(pvalues=pvals[!is.na(pvals)], potential_pvalues=pot_pvals))
}

#' R backend for enrichmentPValue calculation
#' @inheritParams enrichmentPValue
#' @param minorAllele which variant value to count as the minor allele
#' @param famIds family ids corresponding to rows of the snpMap
#' @return p-value
enrichmentPValue_R_Backend <- function(snpMat, famIds, sharingProbs, minorAllele, threshold=0)
{
    # convert matrix to list of families with each allele
    shareList <- convertMatrix(snpMat@.Data, famIds, minorAllele)
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

#' convert snpMatrix to a list of vectors of sharing
#' @inheritParams multipleVariantPValue
#' @param famIds family ids corresponding to rows of the snpMap
#' @return list of boolean vectors indicating sharing pattern for each variant
convertMatrix <- function(snpMat, famIds, minorAllele)
{
    shareList <- sapply(colnames(snpMat), function(var)
    {
        ret <- c()
        for (fid in unique(famIds))
        {
            famRows <- which(famIds == fid)
            famAlleles <- snpMat[famRows, var]
            if (!all(famAlleles != minorAllele[var] & famAlleles != 2))
            {
                ret <- c(ret, all(famAlleles == minorAllele[var] | famAlleles == 2))
                names(ret)[length(ret)] <- fid
            }
        }
        if (length(ret) == 0)
        {
            stop("found a variant not seen in any subject")
        }
        return(ret)
    }, USE.NAMES=TRUE, simplify=FALSE)
    return(shareList[!sapply(shareList, is.null)])
}
