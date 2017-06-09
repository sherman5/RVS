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
    nonFounderNodes <- lapply(setdiff(procPed$id, procPed$founders),
        function(nf) {gRain::cptable(c(nf, procPed$parents[1,nf],
            procPed$parents[2,nf]), values=mendelProbTable, levels=0:2)})

    # create bayesian network
    condProbTable <- gRain::compileCPT(c(founderNodes, nonFounderNodes))
    return(gRain::grain(condProbTable))
}

#' \code{marginalProb} calculates the probability
#'  that a set of nodes are in given states
#'
#' @param net bayesian network from gRain package
#' @param states named list of states for each node
#' @return joint probability
#' @keywords internal
marginalProb <- function(net, states)
{
    prob <- 1
    for (n in names(states))
    {
        if (prob > 0) # prevents conditioning on zero prob events
        {
            # calculate probability for this node
            p <- unname(gRain::querygrain(net, n)[[1]])
            prob <- prob * sum(p[states[[n]]+1])

            # condition on this node being in the correct states
            net <- gRain::setEvidence(net, evidence=sapply(simplify=FALSE,
                X=n, FUN=function(x) as.numeric(0:2 %in% states[[n]])))
        }
    }
    return(prob)
}
