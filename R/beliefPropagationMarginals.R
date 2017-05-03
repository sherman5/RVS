beliefPropagation <- function(graph, fixedVars, marginalVars, prior=NULL)
{
    junctionTree <- jTree(moralize(graph))

    # create all individual pmfs
    vars <- nodes(graph)[!nodes(graph) %in% fixedVars[1,]]
    indivPmf <- list()
    for (v in vars)
    {
        indivPmf[[v]] <- individualPmf(graph, v, fixedVars, prior)
    }

    # place each pmf in a clique containing all variables
    cliques <- sapply(junctionTree$cliques, as.numeric)
    cliquePmfs <- lapply(

    for (c in cliques)
    {
        

    }

    for (v in vars)
    {
        pmf <- indivPmf[[v]]
        valid <- sapply(cliques, function(x) all(pmf[['vars']] %in% x))
        
               
    }

    # mutliply pmfs within cliques

    # pass messages ...

}

largestClique <- function(graph, marginalNodes)
{
    junctionTree <- jTree(moralize(graph))
    return(max(sapply(junctionTree[['cliques']], length)))
}
