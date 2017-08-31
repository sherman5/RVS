# table of probabilties for mendelian inheritance
mendelProbTable <- array(0, c(3,3,3))
x <- expand.grid(p1=0:2, p2=0:2) # parents each having 0,1,2 copies of allele
mendelProbTable[1,,] <- with(x, (2-p1)*(2-p2)/4) # 0 copies in offspring
mendelProbTable[2,,] <- with(x, (p1+p2-p1*p2)/2) # 1 copy in offspring
mendelProbTable[3,,] <- with(x, p1*p2/4)         # 2 copies in offspring

#' create bayesian network from processed pedigree
#' @keywords internal
#'
#' @description Creates a bayesian network using the gRain package.
#' The network is built based on the information in a pedigree object
#' that has been processed using \code{processPedigree}.
#' @param procPed processed Pedigree object
#' @param prior prior on number of alleles for founders
#' @return bayesian network from gRain package
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

#' calculates the marginal probability of a set of nodes
#' @keywords internal
#'
#' @description Given a bayesian network from the gRain package and a 
#' named list of (nodes, states), this function returns the joint-marginal
#' probability of each node taking a value in the specified set of states. 
#' @details This function calculates the probability P(A,B,C) by factoring
#' it into conditional probabilities, i.e. P(A|B,C) * P(B|C) * P(C).
#' Starting at the right side, P(C) is computed and then evidence of C
#' being true is added to the network and P(B) is computed - effectively
#' giving the probability P(B|C). This process continues from right to
#' left until the entire product has been computed.
#' @param net bayesian network from gRain package
#' @param states named list of states for each node
#' @return joint-marginal probability
marginalProb <- function(net, states)
{
    prob <- 1
    for (n in names(states))
    {
        if (prob > 0) # prevents conditioning on zero prob events
        {
            # calculate probability for this node
            p <- unname(gRain::querygrain(net, n, exclude=FALSE)[[1]])
            prob <- prob * sum(p[states[[n]]+1])

            # condition on this node being in the correct states
            net <- gRain::setEvidence(net, evidence=sapply(simplify=FALSE,
                X=n, FUN=function(x) as.numeric(0:2 %in% states[[n]])))
        }
    }
    return(prob)
}
