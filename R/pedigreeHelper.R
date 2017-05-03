validPedigree <- function(ped)
{
    size <- length(ped$id)
    if (size == 0) stop('missing \'id\' field')
    else if (length(ped$findex) != size) stop('findex size invalid')
    else if (length(ped$mindex) != size) stop('mindex size invalid')
    else if (length(ped$affected) != size) stop('affected size invalid')
    else return(TRUE)
}

pedToDAG <- function(ped)
{
    # list used to generate DAG
    edgeList <- list()
    index <- 1

    # process nodes in order
    for (i in 1:length(ped$id))
    {
        edgeList[[index]] <- i
        index <- index + 1
    }

    # find edges in graph
    for (i in 1:length(ped$id))
    {
        offspring <- c(which(ped$findex == i), which(ped$mindex == i))
        for (c in offspring)
        {
            edgeList[[index]] <- c(c, i)
            index <- index + 1
        }
    }
    return (dag(edgeList))
}
