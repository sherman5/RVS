#' \code{oneFounderSharingProb} calculate sharing prob assuming one founder
#'  introduces the variant
#'
#' @param procPed processed Pedigree object
#' @return sharing probability
#' @keywords internal
oneFounderSharingProb <- function(procPed)
{
    # each fraction component of probability
    numer <- denom <- 0

    # set all founders to 0 (no variant)
    net <- createNetwork(procPed)
    net <- gRain::setEvidence(net, as.character(procPed$founders),
        rep('0', length(procPed$founders)))

    # sum over probs, conditioning on each founder introducing variant
    for (f in procPed$founders)
    {
        # condition on founder and calculate distribution
        condNet <- gRain::retractEvidence(net, as.character(f))
        condNet <- gRain::setEvidence(condNet, as.character(f), '1')

        # sum relevant probability       
        numer <- numer + marginalProb(condNet, sapply(simplify=FALSE,
            X=as.character(procPed$carriers), FUN=function(x) 1:2))
        denom <- denom + 1 - marginalProb(condNet, sapply(simplify=FALSE,
            X=as.character(procPed$affected), FUN=function(x) 0))
    }
    return(numer/denom)
}

#' \code{twoFounderSharingProb}
twoFounderSharingProb <- function(procPed)
{
    # set all founders to 0 (no variant)
    net <- createNetwork(procPed)
    net <- gRain::setEvidence(net, as.character(procPed$founders),
        rep('0', length(procPed$founders)))


    weight <- 0.5 #nf * PFU
    numer <- denom <- 0
    remainingFounders <- procPed$founders
    for (f1 in procPed$founders)
    {
        # condition on founder
        net1 <- gRain::retractEvidence(net, as.character(f))
        net1 <- gRain::setEvidence(net1, as.character(f), '1')
       
        for (f2 in remainingFounders)
        {
            # condition on founder
            net2 <- gRain::retractEvidence(net1, as.character(f))
            net2 <- gRain::setEvidence(net2, as.character(f), '1')

            w <- ifelse(f1==f2, weight, 1-weight)

            numer <- numer + w * marginalProb(net2, sapply(simplify=FALSE,
                X=as.character(procPed$carriers), FUN=function(x) 1:2))
            denom <- denom + w * (1 - marginalProb(net2, sapply(simplify=FALSE,
                X=as.character(procPed$affected), FUN=function(x) 0)))

        }
        remainingFounders <- setdiff(remainingFounders, f1)
    }
    return(numer/denom)
}

#' \code{exactSharingProb}
exactSharingProb <- function(procPed, alleleFreq)
{
    prior <- with(data.frame(f=alleleFreq), c((1-f)^2, 2*f*(1-f), f^2))
    net <- createNetwork(procPed, prior)

    numer <- marginalProb(net, sapply(as.character(procPed$carriers),
        function(x) 1:2, simplify=FALSE))
    denom <- 1 - marginalProb(net, sapply(as.character(procPed$affected),
        function(x) 0, simplify=FALSE))
    return (numer/denom)
}

#' \code{RVsharing2} Calculates rare variant sharing probability between
#'  affected subjects in a procPedigree
#'
#' @param ped S3 Pedigree object
#' @param alleleFreq frequency of variant present amount the founders
#' @param kinshipCoeff relationship between founders
#' @param nSimulations number of simulations to run for the monte carlo
#'  approach to calculating sharing probabilities
#' @return sharing probability between all affected subjects in the procPedigree
#' @examples
#'  data("sampleprocPedigrees")
#'  RVsharing2(sampleprocPedigrees[[1]])
#' @export
RVsharing2 <- function(ped, alleleFreq, kinshipCoeff, nSimulations)
{
    # load data for prob calculations
    data(mendelProbTable)

    # pre-process procPedigree
    procPed <- processPedigree(ped)

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
        stop('RVsharing2 does not yet support relatedness among founders')
    }
    else
    {
        return(oneFounderSharingProb(procPed))
    }
}
