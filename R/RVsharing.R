#  [13][10][9]P(12|9,10) [8]P(11|8,13) [7][6]P(6)P(13|6,7)P(10|6,7)
#  [5]P(5) [4]P(9|4,5)P(8|4,5) [3]P(3) [2]P(2)P(7|2,3) [1]P(1)P(4|1,2)

# local
# f(1,2,4)
# f(2,3,7)
# f(3)
# f(4,5,8,9)
# f(5)
# f(6,7,10,13)
# f()
# f(8,11,13)
# f(9,10,12)
# f()
# f()

# global
# f(2,4)
# f(3,4,7)
# f(4,7)
# f(5,7,8,9)
# f(7,8,9)
# f(7,8,9,10,13)
# f(8,9,10,13)
# f(9,10,13)
# f(10,13)
# f(13)
# f()

sumOutVariable <- function(pmf, index)
{
    perm <- 1:length(dim(pmf))
    perm <- perm[perm != index]
    perm <- c(index, perm)
    return(colSums(aperm(pmf, perm)))
}

globalPmf <- function(pmf, vars, vals)
{


}

marginalProbability <- function(ped, marginalNodes, prior)
{
    # get graph info
    graph <- pedToDAG(ped)
    order <- as.numeric(topoSort(graph))
    order <- order[!order %in% marginalNodes[1,]]
    
    #initialize conditional probabilities
    condProb <- function(x,p1,p2) 
    {
        pmf <- c((2-p1)*(2-p2)/4, (p1+p2-p1*p2)/2, p1*p2/4)
        return(pmf[x+1])
    }
    prob <- lapply(1:length(graph@nodes), function(x) condProb)

    # overwrite conditional probs for founders
    for (node in 1:length(order))
    {
        if (sum(getParents(ped, node)) == 0)
        {
            prob[[node]] <- function(x,p1,p2) {prior[[node]][x+1]}
        }
    }

    # create local interaction pmf
    localPmf <- function(nodes, vals)
    {
        prod <- 1
        for (i in nodes)
        {
            prod <- prod * prob[[i]](vals[[i]][1], vals[[i]][2],
                vals[[i]][3])
        }
        return(prod)
    }

    # rolling pmf
    pmf <- list()
    pmf[['vars']] <- c()
    pmf[['prob']] <- array()

    # already processed nodes
    processed <- c()

    # sum over variables in order
    for (node in order)
    {
        # find all un-processed children for this node
        candidates <- c(node, getOffspring(ped, node))
        localVars <- candidates[sapply(candidates, function(x)
            {sum(getParents(ped, x) %in% processed) == 0})]

        # find the parents the localVars depend on
        localVarsDep <- unique(c(sapply(localVars, getParents, ped=ped)))
        localVarsDep <- unlist(localVarsDep[localVarsDep != 0])

        # find the vars required by the new pmf being made (node summed out)
        newVars <- unique(c(localVars, localVarsDep, pmf[['vars']]))
        newVars <- newVars[!newVars %in% marginalNodes[1,]]

        # create new pmf
        newPmf <- list()
        newPmf[['vars']] <- newVars
        newPmf[['prob']] <- array(0, rep(3,length(newVars)))

        # find all variable combinations needed for new pmf
        args <- expand.grid(lapply(1:length(newVars), function(x) 0:2))
        args <- as.matrix(unname(args))
        for (a in 1:nrow(args))
        {
            arg <- args[a,]  
            trios <- list()
            for (v in newVars)
            {
                

            }
            localArg <- arg[1:length(newVars)]
        }

        # sum out node, update pmf, mark node as processed
        newPmf[['prob']] <- sumOutVariable(newPmf[['prob']],
            which(newVars==node))
        pmf <- newPmf
        processed <- c(processed, node)
    }
}

pedToDAG <- function(ped)
{
    edgeList <- list()
    index <- 1

    for (i in 1:getSize(ped))
    {
        edgeList[[index]] <- i
        index <- index + 1
    }

    for (i in 1:getSize(ped))
    {
        for (c in getOffspring(ped, i))
        {
            edgeList[[index]] <- c(c, i)
            index <- index + 1
        }
    }
    return (dag(edgeList))
}

hash <- function(...) 
{
    args <- list(...)
    return (paste(paste(names(args), args, sep='='), collapse='_'))
}

RVsharing <- function(ped, alleleFreq = 0.5)
{
    p <- alleleFreq
    q <- 1 - alleleFreq
    prior <- lapply(1:getSize(ped), function(x) c(p^2, 2*p*q, q^2))

    affected <- which(ped$affected == 1)
    nodes <- matrix(nrow = 2, ncol = length(affected))
    nodes[1,] <- affected
    nodes[2,] <- rep(0, length(affected))
    return (beliefProp(ped, nodes, prior))
}

initPotential <- function(ped, prior)
{
    graph <- pedToDAG(ped)
    potential <- list()
    for (node in as.numeric(nodes(graph)))
    {
        if (sum(getParents(ped, node)) == 0)
        {
            potential[[node]] <- function(x) {return (prior[[node]][x+1])}
        }
        else
        {
            potential[[node]] <- function(x, p1, p2)
            {
                prob <- c((2-p1)*(2-p2)/4, (p1+p2-p1*p2)/2, p1*p2/4)
                return (prob[x+1])
            } 
        }
    }
    return (potential)
}

placeNodesInClique <- function(ped, cliques)
{
    graph <- pedToDAG(ped)
    cliqueNodes <- sapply(1:length(cliques), function(x) NULL)
    for (node in as.numeric(nodes(graph)))
    {
        group <- c(node, getParents(ped, node))
        group <- group[group != 0]
        index <- which(sapply(cliques, function(x) all(group %in% x)))[1]
        cliqueNodes[[index]] <- c(cliqueNodes[[index]], node)
    }
    return (cliqueNodes)
}

beliefProp <- function(ped, nodes, prior)
{
    # create junction tree
    graph <- pedToDAG(ped)
    junctionTree <- jTree(moralize(graph))
    cliques <- lapply(junctionTree[['cliques']], as.numeric)

    ## create individual potentials
    potential <- initPotential(ped, graph, prior)

    ## calculate clique potential
    for (cliq in 1:length(cliquePotential))
    {
        nodes <- cliquePotential[[cliq]]
        varComb <- expand.grid(lapply(1:length(vars), function(x) 0:2))
        varComb <- as.matrix(unname(varComb))
        for (i in 1:nrow(varComb))
        {
            var <- varComb[i,]
            arg <- lapply(nodes, function(x) var[which(nodes == x)])
            print(arg)
        }
        print(cliq)
    }
    return (potential)
}


