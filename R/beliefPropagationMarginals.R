# shafer-shenoy
beliefPropagation <- function(graph, fixedVars, marginalVar, prior=NULL)
{
    # create the junction tree
    juncTree <- jTree(moralize(graph))

    # create all individual pmfs
    vars <- as.numeric(nodes(graph)[!nodes(graph) %in% fixedVars[1,]])
    pmfs <- lapply(vars, function(v) createPmf(graph, v, fixedVars, prior))

    # place each pmf in a clique containing all variables
    cliques <- sapply(juncTree$cliques, as.numeric)
    cliquePmfs <- lapply(cliques, function(c)
        {which(sapply(pmfs, function(p) all(p[['vars']] %in% c)))})
    print(cliquePmfs)
    # mutliply pmfs within cliques
    for (c in 1:length(cliques))
    {
        pmfVars <- cliquePmfs[[c]]
        cliquePmfs[[c]] <- list('prob'=1)
        for (p in pmfVars)
        {
            cliquePmfs[[c]] <- multiplyPmf(cliquePmfs[[c]], pmfs[[p]])
        }
    }

    # pass messages ...
    root <- which(junctionTree$parents == 0)
    treeState <- list('tree'=juncTree, 'pmfs'=cliquePmfs, 'messages'=list())
    collect(root, treeState)
    distribute(root, treeState)
    
    # sum out variables to obtain marginal probability

}

collect <- function(node, treeState)
{
    # recursively collect on all children
    children <- which(junctionTree$parents == node)
    for (c in children)
    {
        treeState <- collect(c, treeState)
    }

    # send message
}

distribute <- function(node, treeState)
{
    parent <- junctionTree$parents[node]
}

largestClique <- function(graph, marginalNodes)
{
    junctionTree <- jTree(moralize(graph))
    return(max(sapply(junctionTree[['cliques']], length)))
}

