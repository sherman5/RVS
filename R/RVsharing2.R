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
    for (f in procPed$founders)
    {
        # condition on founder and calculate distribution
        condNet <- gRain::retractEvidence(net, as.character(f))
        condNet <- gRain::setEvidence(condNet, as.character(f), '1')

        # sum relevant probability       
        numer <- numer + marginalProb(condNet, sapply(simplify=FALSE,
            X=as.character(procPed$carriers), FUN=function(dummy) 1:2))
        denom <- denom + 1 - marginalProb(condNet, sapply(simplify=FALSE,
            X=as.character(procPed$affected), FUN=function(dummy) 0))
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
    cor <- relatedFoundersCorrection(length(procPed$founders), kinshipCoeff)
    numer <- denom <- 0
    remainingFounders <- procPed$founders
    for (f1 in procPed$founders)
    {
        # condition on founder
        net1 <- gRain::retractEvidence(net, as.character(f1))
        net1 <- gRain::setEvidence(net1, as.character(f1), '1')
       
        for (f2 in remainingFounders)
        {
            # condition on founder
            net2 <- gRain::retractEvidence(net1, as.character(f2))
            net2 <- gRain::setEvidence(net2, as.character(f2), '1')

            # weight this prob based on 1 (f1==f2) or 2 (f1!=f2) founders
            w <- ifelse(f1==f2, cor, 1 - cor)

            # calculate conditional probability
            numer <- numer + w * marginalProb(net2, sapply(simplify=FALSE,
                X=as.character(procPed$carriers), FUN=function(dummy) 1:2))
            denom <- denom + w * (1 - marginalProb(net2, sapply(simplify=FALSE,
                X=as.character(procPed$affected), FUN=function(dummy) 0)))

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

    # calculate the conditional probability
    numer <- marginalProb(condNet, sapply(simplify=FALSE,
        X=as.character(procPed$carriers), FUN=function(dummy) 1:2))
    denom <- 1 - marginalProb(condNet, sapply(simplify=FALSE,
        X=as.character(procPed$affected), FUN=function(dummy) 0))
   return (numer/denom)
}

#' \code{RVsharing2} Calculates rare variant sharing probability between
#'  carriers in a pedigree given that the variant in seen in one of the
#   affected subjects
#'
#' @param ped S3 Pedigree object
#' @param alleleFreq allele frequency among the founders
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @param nSimulations number of simulations used in monte carlo calculation
#' @return sharing probability between all carriers in pedigree
#' @examples
#'  data("samplePedigrees")
#'  RVsharing2(samplePedigrees[[1]])
#' @export
RVsharing2 <- function(ped, alleleFreq, kinshipCoeff, nSimulations, carriers)
{
    # load data for prob calculations
    data(mendelProbTable)

    # pre-process procPedigree
    procPed <- processPedigree(ped, carriers)

    # calculate sharing prob with appropiate method
    if (!missing(nSimulations))
    {
        return(monteCarloSharingProb(procPed, alleleFreq, kinshipCoeff,
            nSimulations))
    }
    else if (!missing(alleleFreq))
    {
        return(exactSharingProb(procPed, alleleFreq))
    }        
    else if (!missing(kinshipCoeff))
    {
        return(twoFounderSharingProb(procPed, kinshipCoeff))
    }
    else
    {
        return(oneFounderSharingProb(procPed))
    }
}
