mendelProb <- array(0, c(3,3,3))
args <- expand.grid(p1=0:2, p2=0:2)
mendelProb[1,,] <- with(args, (2-p1) * (2-p2) / 4)
mendelProb[2,,] <- with(args, (p1 + p2 - p1*p2) / 2)
mendelProb[3,,] <- with(args, p1 * p2 / 4)

monteCarloSharingProb <- function(ped, alleleFreq, kinshipCoeff,
nSimulations)
{
    denom <- numer <- 0
    defStates <- rep(NA, length(ped$ped$id))
    nFounders <- length(ped$founders)

    if (!missing(alleleFreq))
    {
        p <- with(data.frame(f=alleleFreq), c((1-f)^2, 2*f*(1-f), f^2))
        for (n in 1:nSimulations)
        {
            states <- defStates
            states[ped$founders] <- sample.int(3,nFounders,TRUE,p)-1
            res <- simulateTree(ped, states)

            if (sum(res) >= 1) denom <- denom + 1
            if (all(res >= 1)) numer <- numer + 1
        }
    }
    else if (!missing(kinshipCoeff))
    {
    
    }
    else # one founder introduces
    {
        for (n in 1:nSimulations)
        {
            states <- defStates
            states[ped$founders] <- 0
            states[sample(ped$founders,1)] <- 1
            res <- simulateTree(ped, states)

            if (sum(res) >= 1) denom <- denom + 1
            if (all(res >= 1)) numer <- numer + 1
        }
    }
    return(numer/denom)
}

simulateTree <- function(ped, states)
{
    remain <- which(is.na(states))
    while (length(remain))
    {
        for (i in remain)
        {
            p1 <- states[ped$parents[,i][1]]
            p2 <- states[ped$parents[,i][2]]
        
            if (!is.na(p1) & !is.na(p2))
            {
                states[i] <- sample.int(3,1,p=mendelProb[,p1+1,p2+1])-1
            }
        }
        remain <- which(is.na(states))
    }
   return(states[ped$affected])
}
