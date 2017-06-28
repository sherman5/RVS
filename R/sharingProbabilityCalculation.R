numerProb <- function(net, procPed)
{
    rvInCarriers <- sapply(simplify=FALSE, FUN=function(dummy) 1:2,
        X=as.character(procPed$carriers))
    noRvInNonCarriers <- sapply(simplify=FALSE, FUN=function(dummy) 0,
        X=as.character(setdiff(procPed$affected, procPed$carriers)))
    return(marginalProb(net, c(rvInCarriers, noRvInNonCarriers)))
}

denomProb <- function(net, procPed)
{
    noRvInAny <- sapply(simplify=FALSE, FUN=function(dummy) 0,
        X=as.character(procPed$affected))
    return(1 - marginalProb(net, noRvInAny))
}

#' \code{oneFounderSharingProb} calculate sharing probability assuming one
#'  founder introduces the variant
#'
#' @param procPed pedigree that has been through processPedigree()
#' @return sharing probability
#' @keywords internal
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

#' \code{twoFounderSharingProb} calculate the sharing probability assuming
#'  that at most two founders introduce the variant, accounting for 
#'  relatedness among founders
#'
#' @param procPed pedigree that has been through processPedigree()
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @return sharing probability
#' @keywords internal
twoFounderSharingProb <- function(procPed, kinshipCoeff)
{
    # set all founders to 0 (no variant)
    net <- createNetwork(procPed)
    net <- gRain::setEvidence(net, as.character(procPed$founders),
        rep('0', length(procPed$founders)))

    # calculate kinship correction and initialize variables
    numer <- denom <- 0
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
            w <- ifelse(f1==f2, cor, 1 - cor)
            numer <- numer + w * numerProb(net2, procPed)
            denom <- denom + w * denomProb(net2, procPed)
        }
        remainingFounders <- setdiff(remainingFounders, f1)
    }
    return(numer/denom)
}

#' \code{exactSharingProb} calculates the exact sharing probability given
#'  allele frequency among the founders
#'
#' @param procPed pedigree that has been through processPedigree()
#' @param alleleFreq allele frequency among the founders
#' @return sharing probability
#' @keywords internal
exactSharingProb <- function(procPed, alleleFreq)
{
    # create network with founders having a prior based on the allele freq
    prior <- with(data.frame(f=alleleFreq), c((1-f)^2, 2*f*(1-f), f^2))
    net <- createNetwork(procPed, prior)

    # compute probability
    return (numerProb(net, procPed) / denomProb(net, procPed))
}
