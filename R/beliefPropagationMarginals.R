beliefPropagationMarginals <- function(graph, marginalNodes, prior)
{
    junctionTree <- jTree(moralize(graph))

    # create all individual pmfs

    # place each pmf in a clique containing all variables

    # mutliply pmfs withing cliques

    # pass messages ...


}

largestClique <- function(graph, marginalNodes)
{
    junctionTree <- jTree(moralize(graph))
    return(max(sapply(junctionTree[['cliques']], length)))
}
