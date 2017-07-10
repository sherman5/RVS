#' numerator of sharing probability
#' @keywords internal
#'
#' @description calculates the numerator of the sharing probability
#'  outline in section 2.1 of Bureau et al.
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
#'  outline in section 2.1 of Bureau et al.
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
#'  the pedigree. Condition on each founder and sum over all resulting
#'  probabilities. 
#' @param procPed pedigree that has been through processPedigree()
#' @return sharing probability
oneFounderSharingProb <- function(procPed)
{
    # set all founders to 0 (no variant)
    net <- createNetwork(procPed)
    net <- gRain::setEvidence(net, as.character(procPed$founders),
        rep('0', length(procPed$founders)))

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

#' calculate the sharing probability when two founders can introduce variant
#' @keywords internal
#'
#' @description In the case of relatedness among founders, assume that up
#'  to two founders could introduce the variant and condition on all possible
#'  pairs.
#' @param procPed pedigree that has been through processPedigree()
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @return sharing probability
twoFounderSharingProb <- function(procPed, kinshipCoeff)
{
    # set all founders to 0 (no variant)
    net <- createNetwork(procPed)
    net <- gRain::setEvidence(net, as.character(procPed$founders),
        rep('0', length(procPed$founders)))

    # calculate kinship correction and initialize variables
    numer1 <- denom1 <- numer2 <- denom2 <- 0
    remainingFounders <- procPed$founders
    cor <- relatedFoundersCorrection(length(procPed$founders), kinshipCoeff)
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
#            w <- ifelse(f1==f2, cor, 1 - cor)
#            numer <- numer + w * numerProb(net2, procPed)
#            denom <- denom + w * denomProb(net2, procPed)
            if (f1 == f2)
            {
                numer1 <- numer1 + numerProb(net2, procPed)
                denom1 <- denom1 + denomProb(net2, procPed)
            }
            else
            {
                numer2 <- numer2 + numerProb(net2, procPed)
                denom2 <- denom2 + denomProb(net2, procPed)
            }
        }
        remainingFounders <- setdiff(remainingFounders, f1)
    }
    return(cor * numer1/denom1 + (1 - cor) * numer2/denom2)
}

#' exact sharing probability calculation
#' @keywords internal
#'
#' @description Calculate the exact sharing probability given the minor allele
#'  frequency among the founders (population).
#' @param procPed pedigree that has been through processPedigree()
#' @param alleleFreq allele frequency among the founders
#' @return sharing probability
exactSharingProb <- function(procPed, alleleFreq)
{
    # create network with founders having a prior based on the allele freq
    prior <- with(data.frame(f=alleleFreq), c((1-f)^2, 2*f*(1-f), f^2))
    net <- createNetwork(procPed, prior)

    # compute probability
    return (numerProb(net, procPed) / denomProb(net, procPed))
}
