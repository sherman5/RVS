#' \code{monteCarloSharingProb}
monteCarloSharingProb <- function(procPed, alleleFreq, kinshipCoeff,
nSimulations)
{
    if (!missing(alleleFreq))
    {
        p <- with(data.frame(f=alleleFreq), c((1-f)^2, 2*f*(1-f), f^2))
        founderFunc <- function(n) sample.int(3,n,TRUE,p) - 1
    }
    else if (!missing(kinshipCoeff))
    {
        stop('RVsharing2 does not yet support relatedness among founders')
    }
    else # one founder introduces
    {
        founderFunc <- function(n) sample(c(rep(0,n-1),1))
    }
    return(runMonteCarlo(procPed, founderFunc, nSimulations))
}

#' \code{runMonteCarlo}
runMonteCarlo <- function(procPed, founderFunc, nSim)
{
    numer <- denom <- 0
    defaultStates <- rep(NA, procPed$size)

    for (n in 1:nSim) #TODO: run simulations in parallel
    {
        states <- defaultStates
        states[procPed$founders] <- founderFunc(length(procPed$founders))   
        res <- simulatePedigree(procPed, states)
        
        if (sum(res) >= 1) denom <- denom + 1
        if (all(res >= 1)) numer <- numer + 1
    }
    return(numer/denom)
}

#' \code{simulatePedigree}
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
