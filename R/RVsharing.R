sumOutVariable <- function(pmf, index)
{
    perm <- 1:length(dim(pmf))
    perm <- perm[perm != index]
    perm <- c(index, perm)
    return(colSums(aperm(pmf, perm)))
}

marginalProbability <- function(ped, nodes, prior)
{
    # get graph info
    graph <- pedToDAG(ped)
    order <- as.numeric(topoSort(graph))
    print(order)
    edges <- edgeMatrix(graph)
    
    #initialize conditional probabilities
    condProb <- function(x,p1,p2) 
    {
        pmf <- c((2-p1)*(2-p2)/4, (p1+p2-p1*p2)/2, p1*p2/4)
        return(pmf[x+1])
    }
    prob <- lapply(1:length(order), function(x) condProb)

    # overwrite conditional probs for founders
    for (node in 1:length(order))
    {
        if (sum(getParents(ped, node)) == 0)
        {
            prob[[node]] <- function(x,p1,p2) {prior[[node]][x+1]}
        }
    }

    # create local interaction pmf
    localPmf <- function(nodes, args)
    {
        prod <- 1
        for (i in 1:length(nodes))
        {
            prod <- prod * prob[[nodes[i]]](args[[i]][1], args[[i]][2],
                args[[i]][3])
        }
        return(prod)
    }

    # sum over variables in order
    vars <- c()
    pmf <- array()
    processed <- c()
    for (node in order)
    {
        # find all interactions (children & spouses) for this node
        children <- sapply(getOffspring(ped, node), function(x)
            {sum(getParents(ped, x) %in% processed) == 0})
        parents <- sapply(children, getParents, ped = ped)

        # include this node if a parent was not previously processed
        allNodes <- children
        if (sum(getParents(ped, node) %in% processed) == 0)
            {allNodes <- c(allNodes, node)}

        # create current pmf (function)
        mainVars <- as.numeric(unlist(unique(c)))
        allVars <- as.numeric(unlist(unique(c(allNodes, c(parents)))))

        # record primary vars and parents, set node to 0,1,2
        # sum over 3 values of node to get values for pmf

        # update main pmf
        newPmf <- array(0, rep(3, length(allVars)))
        args <- expand.grid(lapply(1:length(allVars), function(x) 0:2))
        args <- as.matrix(unname(args)) + 1
        for (a in 1:nrow(args))
        {
            arg <- args[a,]
            curArg <- arg[1:length(curVars)]
            oldArg <- arg[(length(curVars)+1):length(oldVars)]
#            arg <- matrix(args[a,], 1)
            newPmf[arg] <- localPmf(curArg) * pmf[matrix(oldArg, 1)]
        }
l
        # sum out node
        pmf <- sumOutVariable(newPmf, index(node))

        # mark this node as processed
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


