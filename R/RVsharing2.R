#' \code{createNetwork} create bayesian network from procPedigree
#'
#' @param procPed processed Pedigree object
#' @param prior prior on number of alleles for founders
#' @return bayesian network friom gRain package
#' @keywords internal
createNetwork <- function(procPed, prior=c(1,2,1))
{
    # process founders
    founderNodes <- lapply(procPed$founders, gRain::cptable, values=prior,
        levels=0:2)

    # process non-founders
    f <- function(nf)
    {
        gRain::cptable(c(nf, procPed$parents[1,nf], procPed$parents[2,nf]),
            values=mendelProbTable, levels=0:2)
    }
    nonFounderNodes <- lapply(setdiff(1:procPed$size, procPed$founders), f)

    # create bayesian network
    condProbTable <- gRain::compileCPT(c(founderNodes, nonFounderNodes))
    return(gRain::grain(condProbTable))
}

#' \code{marginalProb} calculates probability that 1) all marginal nodes
#'  are 0 and 2) all marginal nodes are 1
#'
#' @param net Bayesian network from gRain package
#' @param marginalNodes nodes in the joint-marginal distribution
#' @return two probabilities
#' @keywords internal
marginalProb <- function(net, marginalNodes)
{
    p0 <- p1 <- 1
    net0 <- net1 <- net
    for (n in as.character(marginalNodes))
    {
        if (p0 > 0) # prevents conditioning on zero prob events
        {
            p0 <- p0 * unname(gRain::querygrain(net0, n)[[1]][1])
            net0 <- gRain::setEvidence(net0, n, '0')
        }
        if (p1 > 0)
        {
            p1 <- p1 * unname(gRain::querygrain(net1, n)[[1]][2])
            net1 <- gRain::setEvidence(net1, n, '1')
        }
    }
    return(c(p0,p1))
}

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
        net <- gRain::retractEvidence(net, as.character(f))
        net <- gRain::setEvidence(net, as.character(f), '1')
        prob <- marginalProb(net, procPed$affected)

        # sum relevant probability       
        numer <- numer + prob[2]
        denom <- denom + 1 - prob[1]

        # reset founder
        net <- gRain::retractEvidence(net, as.character(f))
        net <- gRain::setEvidence(net, as.character(f), '0')
    }
    return(numer/denom)
}

#' \code{exactSharingProb}
exactSharingProb <- function(procPed, alleleFreq)
{
    p <- with(data.frame(f=alleleFreq), c((1-f)^2, 2*f*(1-f), f^2))
    net <- createNetwork(pprocPed, p)
    prob <- marginalProb(net, procPed$affected)
    return (prob[2] / (1 - prob[1]))
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
    # pre-process procPedigree
    procPed <- processprocPedigree(procPed)

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
