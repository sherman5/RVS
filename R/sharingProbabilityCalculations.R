#' numerator of sharing probability
#' @keywords internal
#'
#' @description calculates the numerator of the sharing probability
#' outline in section 2.1 of Bureau et al.
#' @param gRain bayesian network
#' @param procPed pedigree object that has been process with processPedigree
#' @return numerator value
numerProb <- function(net, procPed)
{
    rvInCarriers <- sapply(simplify=FALSE, FUN=function(dummy) 1:2,
        X=as.character(procPed$carriers))
    noRvInNonCarriers <- sapply(simplify=FALSE, FUN=function(dummy) 0,
        X=as.character(setdiff(procPed$affected, procPed$carriers)))
    return(marginalProb(net, c(rvInCarriers, noRvInNonCarriers)))
}

#' denominator of sharing probability
#' @keywords internal
#'
#' @description calculates the denominator of the sharing probability
#' outline in section 2.1 of Bureau et al.
#' @param gRain bayesian network
#' @param procPed pedigree object that has been process with processPedigree
#' @return denominator value
denomProb <- function(net, procPed)
{
    noRvInAny <- sapply(simplify=FALSE, FUN=function(dummy) 0,
        X=as.character(procPed$affected))
    return(1 - marginalProb(net, noRvInAny))
}

#' calculate sharing probability in basic case
#' @keywords internal
#'
#' @description Assume that only one founder can introduce the variant to 
#' the pedigree. Condition on each founder and sum over all resulting
#' probabilities. 
#' @param procPed pedigree that has been through processPedigree()
#' @return sharing probability
oneFounderSharingProb <- function(procPed)
{
    # set all founders to 0 (no variant)
    net <- try(createNetwork(procPed))
    if (class(net)[1]=="try-error") stop("Creation of Bayesian network for the pedigree failed.")
    net <- try(gRain::setEvidence(net, as.character(procPed$founders),
        rep('0', length(procPed$founders))))
    if (class(net)[1]=="try-error") stop("Setting founder genotypes of the pedigree failed.")
    
    # sum over probs, conditioning on each founder introducing variant
    numer <- denom <- 0
    for (f in procPed$founders) #TODO: use sapply here
    {
        # condition on founder and calculate distribution
        condNet <- gRain::retractEvidence(net, as.character(f))
        condNet <- gRain::setEvidence(condNet, as.character(f), '1')

        # compute probability
        denom <- denom + denomProb(condNet, procPed)
        numer <- numer + numerProb(condNet, procPed)
    }
    return(numer/denom)
}

#' sharing probability when founder pair introduces variant
#' @keywords internal
#'
#' @description In the case of relatedness among founders, assume that up
#' to two founders could introduce the variant and condition on all possible
#' pairs.
#' @param procPed pedigree that has been through processPedigree()
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @param kinshipOrder order of the polynomial approximation to the distribtion
#' of the number of distinct alleles in the founders (d in Bureau et al.).
#' Must be <= 5
#' @return sharing probability
twoFounderSharingProb <- function(procPed, kinshipCoeff, kinshipOrder)
{
    # set all founders to 0 (no variant)
    net <- try(createNetwork(procPed))
    if (class(net)[1]=="try-error") stop("Creation of Bayesian network for the pedigree failed.")
    net <- try(gRain::setEvidence(net, as.character(procPed$founders),
        rep('0', length(procPed$founders))))
    if (class(net)[1]=="try-error") stop("Setting founder genotypes of the pedigree failed.")
    
    # calculate kinship correction and initialize variables
    numer <- denom <- 0
    remainingFounders <- procPed$founders
    cor <- relatedFoundersCorrection(length(procPed$founders),
        kinshipCoeff, kinshipOrder)
    for (f1 in procPed$founders) #TODO: use sapply here
    {
        # condition on founder
        net1 <- gRain::retractEvidence(net, as.character(f1))
        net1 <- gRain::setEvidence(net1, as.character(f1), '1')

        for (f2 in remainingFounders)
        {
            # condition on founder
            net2 <- gRain::retractEvidence(net1, as.character(f2))
            net2 <- gRain::setEvidence(net2, as.character(f2), '1')

            # calculate (weighted) conditional probability
            # Correction to formula 2 of Bioinformatics paper
            w <- ifelse(f1==f2, cor, 2*(1-cor) / (length(procPed$founders)-1))
            numer <- numer + w * numerProb(net2, procPed)
            denom <- denom + w * denomProb(net2, procPed)
        }
        remainingFounders <- setdiff(remainingFounders, f1)
    }
    return(numer/denom)
}

#' exact sharing probability calculation
#' @keywords internal
#'
#' @description Calculate the exact sharing probability given the minor allele
#' frequency among the founders (population).
#' @param procPed pedigree that has been through processPedigree()
#' @param alleleFreq allele frequency among the founders
#' @return sharing probability
exactSharingProb <- function(procPed, alleleFreq)
{
    # create network with founders having a prior based on the allele freq
    prior <- c((1-alleleFreq)^2, 2*alleleFreq*(1-alleleFreq), alleleFreq^2)
    net <- createNetwork(procPed, prior)

    # compute probability
    return (numerProb(net, procPed) / denomProb(net, procPed))
}
