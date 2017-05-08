RVsharing <- function(ped)
{
    affected <- which(ped$affected == 1)
    margNodes <- matrix(nrow=2, ncol=length(affected))
    margNodes[1,] <- affected
    prior <- lapply(1:length(fam$id), function(x) c(1,0,0))
        
    graph <- pedToDAG(ped)
    totalParents <- sapply(1:length(nodes(graph)),
        function(x) sum(as.numeric(parents(x, graph))))
    founders <- which(totalParents == 0)

    numer <- 0
    denom <- 0
    for (f in founders)
    {
        prior[[f]] <- c(0,1,0)
        margNodes[2,] <- rep(0, length(affected))
        margProb_0 <- variableEliminationMarginals(graph, margNodes, prior)
        margNodes[2,] <- rep(1, length(affected))
        margProb_1 <- variableEliminationMarginals(graph, margNodes, prior)
        prior[[f]] <- c(1,0,0)

        numer <- numer + margProb_1
        denom <- denom + 1 - margProb_0
    }
    return(numer/denom)
}

