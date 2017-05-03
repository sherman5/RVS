variableEliminationMarginals <- function(graph, marginalNodes, prior)
{
    # get graph order, max size of pmf
    order <- as.numeric(topoSort(graph))
    order <- order[!order %in% marginalNodes[1,]]
    maxSize <- largestPmf(graph, marginalNodes)

    # partial sum pmf
    pmf <- list()
    pmf[['vars']] <- c()
    pmf[['prob']] <- 1

    # already processed nodes
    processed <- c()

    # sum over variables in order
    for (node in order)
    {
        # find all un-processed children for this node
        candidates <- c(node, as.numeric(children(node,graph)))
        vars <- candidates[sapply(candidates, function(x)
            {sum(as.numeric(parents(x, graph)) %in% processed) == 0})]

        # find the parents that vars depend on, create local pmf
        localPmf <- createLocalPmf(graph, vars, prior, marginalNodes)

        # update old pmf with new local pmf
        pmf <- multiplyPmf(pmf, localPmf)

        # sum out node, update pmf, mark node as processed
        pmf <- sumOutVariable(pmf, node)
        processed <- c(processed, node)
    }
    return (pmf[['prob']][1])
}

largestPmf <- function(graph, marginalNodes)
{
    # create order
    order <- as.numeric(topoSort(graph))
    order <- order[!order %in% marginalNodes[1,]]

    # check the pmf associated with each node
    processed <- c()
    current <- c()
    maxSize <- 0
    for (node in order)
    {
        # find all un-processed children for this node
        candidates <- c(node, as.numeric(children(node,graph)))
        vars <- candidates[sapply(candidates, function(x)
            {sum(as.numeric(parents(x, graph)) %in% processed) == 0})]

        # find dependencies of the current variables
        varsDep <- sapply(vars, function(x) as.numeric(parents(x,graph)))
        current <- unique(c(current, vars, unlist(c(varsDep))))
        current <- current[current %in% order]
        maxSize <- max(maxSize, length(current))

        # mark node as processed, remove for rolling list of vars
        processed <- c(processed, node)
        current <- current[current != node]
    }
    return(maxSize)   
}
