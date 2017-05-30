#' \code{createNetwork} create bayesian network from pedigree
#'
#' @param ped S3 Pedigree object
#' @param prior prior on number of alleles for founders
#' @return bayesian network friom gRain package
#' @keywords internal
createNetwork <- function(ped, prior=c(1,1,1))
{
    # process founders
    founderNodes <- lapply(ped$founders, gRain::cptable, values=prior,
        levels=0:2)

    # process non-founders
    f <- function(nf)
    {
        gRain::cptable(c(nf, ped$parents[1,nf], ped$parents[2,nf]),
            values=mendelProb, levels=0:2)
    }
    nonFounderNodes <- lapply(setdiff(ped$ped$id, ped$founders), f)

    # create bayesian network
    return(gRain::grain(gRain::compileCPT(c(founderNodes,nonFounderNodes))))
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
#' @param ped pedigree S3 Object
#' @return sharing probability
#' @keywords internal
oneFounderSharingProb <- function(ped)
{
    # each fraction component of probability
    numer <- denom <- 0

    # set all founders to 0 (no variant)
    net <- createNetwork(ped)
    net <- gRain::setEvidence(net, as.character(ped$founders),
        rep('0', length(ped$founders)))

    # sum over probs, conditioning on each founder introducing variant
    for (f in ped$founders)
    {
        # condition on founder and calculate distribution
        net <- gRain::retractEvidence(net, as.character(f))
        net <- gRain::setEvidence(net, as.character(f), '1')
        prob <- marginalProb(net, ped$affected)

        # sum relevant probability       
        numer <- numer + prob[2]
        denom <- denom + 1 - prob[1]

        # reset founder
        net <- gRain::retractEvidence(net, as.character(f))
        net <- gRain::setEvidence(net, as.character(f), '0')
    }
    return(numer/denom)
}

exactSharingProb <- function(ped, alleleFreq)
{
}

#' \code{RVsharing2} Calculates rare variant sharing probability between
#'  affected subjects in a pedigree
#'
#' @param ped pedigree object (S3)
#' @param alleleFreq frequency of variant present amount the founders
#' @param kinshipCoeff relationship between founders
#' @param nSimulations number of simulations to run for the monte carlo
#'  approach to calculating sharing probabilities
#' @return sharing probability between all affected subjects in the pedigree
#' @examples
#'  data("samplePedigrees")
#'  RVsharing2(samplePedigrees[[1]])
#' @export
RVsharing2 <- function(ped, alleleFreq, kinshipCoeff, nSimulations)
{
    # pre-process pedigree
    ped <- processPedigree(ped)

    # calculate sharing prob with appropiate method
    if (!missing(nSimulations))
    {
        return(monteCarloSharingProb(ped, alleleFreq, kinshipCoeff,
            nSimulations))
    }
    else if (!missing(alleleFreq))
    {
        return(exactSharingProb(ped, alleleFreq))
    }        
    else if (!missing(kinshipCoeff))
    {
        stop('RVsharing2 does not yet support relatedness among founders')
    }
    else
    {
        return(oneFounderSharingProb(ped))
    }
}
