#' @include grainNetworkHelper.R
NULL

#' calculates sharing probability by simulating pedigree outcomes
#' @keywords internal
#'
#' @description Calculates the same exact probability as RVsharing,
#'  except uses monte carlo simulation instead of exact computation.
#'  This method allows for more flexibility in the scenarios considered.
#' @param procPed pedigree that has been through \code{processPedigree}
#' @inheritParams RVsharing
#' @return sharing probability between all carriers in pedigree
monteCarloSharingProb <- function(procPed, alleleFreq, kinshipCoeff,
nSim, founderDist)
{
    if (!missing(alleleFreq)) # known allele frequency in population
    {
        p <- with(data.frame(f=alleleFreq), c((1-f)^2, 2*f*(1-f), f^2))
        founderDist <- function(n) sample.int(3,n,TRUE,p) - 1
    }
    else if (!missing(kinshipCoeff)) # related founders
    {
        w <- relatedFoundersCorrection(length(procPed$founders),
            kinshipCoeff) # prob one founder introduces variant
        founderDist <- function(n)
            {sample(c(rep(0,n-2), 1, ifelse(runif(1) < w, 0, 1)))}
    }
    else if (missing(founderDist)) # one founder introduces
    {
        founderDist <- function(n) sample(c(rep(0,n-1),1))
    }
    return(runMonteCarlo(procPed, founderDist, nSim))
}

#' run the monte carlo simulation
#' @keywords internal
#'
#' @description Given a number of simulations and a distribution
#'  of variants in the founders, this function simulates possbile
#'  outcomes of the pedigree and returns a sharing probability.
#' @details If the number of simulations is greater than 20,000
#'  then the computation is done in parallel (as long as the package
#'  parallel is available)
#' @inheritParams monteCarloSharingProb
#' @return sharing probability between all carriers in pedigree
runMonteCarlo <- function(procPed, founderDist, nSim)
{
    oneSim <- function(dummy) # carry out one simulation of the pedigree
    {
        states <- rep(NA, procPed$size)
        states[procPed$founders] <- founderDist(length(procPed$founders))   
        sim <- simulatePedigree(procPed, states)
        return(c(all(sim$carriers >= 1), sum(sim$affected) >= 1))
    }

    if (is.element('parallel', installed.packages()[,1]) & nSim > 2e4)
    {
        nCores <- parallel::detectCores()
        print(paste('RVsharing running in parallel with:', nCores, 'cores'))
        cl <- parallel::makeCluster(nCores)
        prob <- parallel::parSapply(cl, 1:nSim, oneSim)
        parallel::stopCluster(cl)
    }
    else
    {
        prob <- sapply(1:nSim, oneSim)
    }
    variantPresent <- prob[2,]
    return(sum(prob[1,variantPresent]) / sum(variantPresent))
}

#' simulates pedigree given founder states
#' @keywords internal
#'
#' @description Given the states (number of allele copies) of the founders,
#'  this function simulates mendelian inheritance and returns the states
#'  of all subjects in the pedigree
#' @param procPed pedigree that has been through processPedigree()
#' @param states state of each founder (0,1,2 copies of variant)
#' @return states for all subjects in pedigree
simulatePedigree <- function(procPed, states)
{
    remain <- which(is.na(states)) # undecided states
    while (length(remain)) # loop until all states decided
    {
        for (i in remain) # check if each state is ready to be decided
        {
            p1 <- states[procPed$parents[1,i]] # find parent states
            p2 <- states[procPed$parents[2,i]]

            if (!is.na(p1) & !is.na(p2)) # both parents decided
            {
                states[i] <- sample.int(3, 1,
                    prob=mendelProbTable[,p1+1,p2+1]) - 1 # decide state
            }
        }
        remain <- which(is.na(states)) # update which states remaining
    }
    return(list(affected=states[procPed$affected],
        carriers=states[procPed$carriers]))
}

#' depreciated function
#' @export
#' @rdname GeneDrop
#'
#' @description This function is depreciated with version >= 2.0
#'  and should not be used, instead use RVsharing with nSim option
#' @param ... arguments to the old function
#' @return none
GeneDrop <- function(...)
{
    stop(paste('function depreciated with version >= 2.0, use the',
        '\'RVsharing\' function with option \'nSim\''))
}

#' @export
#' @rdname GeneDrop
GeneDropSim.allsubsets.fn <- GeneDrop

#' @export
#' @rdname GeneDrop
GeneDropSim.fn <- GeneDrop

#' @export
#' @rdname GeneDrop
GeneDropSimExcessSharing.fn <- GeneDrop


