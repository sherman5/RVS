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
    oneSim <- function(x)
    {
        states <- rep(NA, procPed$size)
        states[procPed$founders] <- founderFunc(length(procPed$founders))   
        res <- simulatePedigree(procPed, states)
        return(c(all(res >= 1), sum(res) >= 1))
    }

    if (is.element('parallel', installed.packages()[,1]) & nSim > 2e4)
    {
        cl <- parallel::makeCluster(parallel::detectCores())
        parallel::clusterExport(cl, 'mendelProbTable')
        prob <- parallel::parSapply(cl, 1:nSim, oneSim)
        parallel::stopCluster(cl)
    }
    else
    {
        prob <- sapply(1:nSim, oneSim)
    }
    return(sum(prob[1,]) / sum(prob[2,]))
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
