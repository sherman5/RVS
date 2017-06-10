#' \code{monteCarloSharingProb} calculates sharing probability by
#'  simulating many outcomes of the pedigree
#'
#' @param procPed pedigree that has been through processPedigree()
#' @param alleleFreq allele frequency among the founders
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @param nSimulations number of simulations used in calculation
#' @return sharing probability
#' @keywords internal
monteCarloSharingProb <- function(procPed, founderDist, alleleFreq, kinshipCoeff, nSimulations)
{
    if (!missing(alleleFreq))
    {
        p <- with(data.frame(f=alleleFreq), c((1-f)^2, 2*f*(1-f), f^2))
        founderDist <- function(n) sample.int(3,n,TRUE,p) - 1
    }
    else if (!missing(kinshipCoeff))
    {
        w <- relatedFoundersCorrection(length(procPed$founders),
            kinshipCoeff)
        founderDist <- function(n)
        {
            if (runif(1) < w)
                sample(c(rep(0,n-1),1)) # one founder introduces
            else
                sample(c(rep(0,n-2),c(1,1))) # two founders introduce
        }        
    }
    else if (missing(founderDist))# one founder introduces
    {
        founderDist <- function(n) sample(c(rep(0,n-1),1))
    }
    return(runMonteCarlo(procPed, founderDist, nSimulations))
}

#' \code{runMonteCarlo} run the monte carlo simulation
#'
#' @param procPed pedigree that has been through processPedigree()
#' @param founderFunc function describing how to sample alleles
#'  for the founders
#' @param nSimulations number of simulations used in calculation
#' @return sharing probability
#' @keywords internal
runMonteCarlo <- function(procPed, founderDist, nSimulations)
{
    oneSim <- function(x)
    {
        states <- rep(NA, procPed$size)
        states[procPed$founders] <- founderDist(length(procPed$founders))   
        res <- simulatePedigree(procPed, states)
        return(c(all(res >= 1), sum(res) >= 1))
    }

    if (is.element('parallel', installed.packages()[,1]) & nSimulations>2e4)
    {
        cl <- parallel::makeCluster(parallel::detectCores())
        parallel::clusterExport(cl, 'mendelProbTable')
        prob <- parallel::parSapply(cl, 1:nSimulations, oneSim)
        parallel::stopCluster(cl)
    }
    else
    {
        prob <- sapply(1:nSimulations, oneSim)
    }
    return(sum(prob[1,]) / sum(prob[2,]))
}

#' \code{simulatePedigree} simulates pedigree given founder states
#'
#' @param procPed pedigree that has been through processPedigree()
#' @param states state of each founder (0,1,2 copies of variant)
#' @return states for all subjects in pedigree
#' @keywords internal
simulatePedigree <- function(procPed, states)
{
    remain <- which(is.na(states))
    while (length(remain))
    {
        for (i in remain)
        {
            p1 <- states[procPed$parents[1,i]]
            p2 <- states[procPed$parents[2,i]]

            if (!is.na(p1) & !is.na(p2))
            {
                states[i] <- sample.int(3, 1,
                    prob=mendelProbTable[,p1+1,p2+1]) - 1
            }
        }
        remain <- which(is.na(states))
    }
    return(states[procPed$affected])
}
