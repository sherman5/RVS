getSize <- function(ped)
{
    return (length(ped$id))
}

getOffspring <- function(ped, index)
{
    return (c(which(ped[['findex']] == index),
        which(ped[['mindex']] == index)))    
}

getParents <- function(ped, index)
{
    return (c(ped$findex[index], ped$mindex[index]))
}

getSpouses <- function(ped, index)
{
    # get all kids
    kids <- getOffspring(ped, index)
    if (length(kids) == 0) {return (0)}

    # find all parents
    parents <- sapply(kids, getParents, ped = ped)

    # remove index from list of parents, keep unique indices
    return (unique(parents[which(parents != index)]))
}

getInLaws <- function(ped, index)
{
    spouses <- getSpouses(ped, index)
    if (spouses[1] == 0) {return (c(0))} # unmarried

    inLaws <- sapply(spouses, getParents, ped = ped)
    inLaws <- unique(c(inLaws))
    if (length(inLaws) == 1) {inLaws <- c(0,0)} # married to founder
    return (inLaws)
}

getFounders <- function(ped)
{
    sumParents <- function(i)
    {
        return (sum(getParents(ped, i)))
    }

    totalParents <- sapply(1:getSize(ped), sumParents)
    return (which(totalParents == 0))
}

getFinalDescendants <- function(ped)
{
    sumOffspring <- function(i) {return (sum(getOffspring(ped, i)))}
    totalKids <- sapply(1:getSize(ped), sumOffspring)
    return (which(totalKids == 0))
}

getNextGeneration <- function(ped, generation)
{
    allChildren <- sapply(generation, getOffspring, ped = ped)
    return (unique(c(unlist(allChildren))))
}

getDescendants <- function(ped, founder)
{
    # used to label descendant depth
    labelDepth <- function(x, d) {return (c(x, d))}

    # get first generation
    gen <- getNextGeneration(ped, founder)
    if (length(gen) == 0) {return (matrix(nrow = 0, ncol = 2))}

    # process each generation
    depth <- 1
    descend <- matrix(nrow = 2, ncol = 0)
    while (length(gen) > 0)
    {
        descend <- cbind(descend, sapply(gen, labelDepth, d = depth))
        gen <- getNextGeneration(ped, gen)
        depth <- depth + 1
    }
    return (descend)
}

getTopFounders <- function(ped)
{
    sumParents <- function(i)
    {
        return (sum(getParents(ped, i) + sum(getInLaws(ped, i))))
    }
    totalParents <- sapply(1:getSize(ped), sumParents)
    return (which(totalParents == 0))
}

# depth first search of parents going back to founder 
getLineage <- function(ped, founder, descendant)
{
    parents <- getParents(ped, descendant)
    allDes <- getDescendants(ped, founder)[1,]    
    
    if (!(founder %in% parents) && descendant %in% allDes)
    {
        lin <- parents
        if (parents[1] %in% allDes)
        {
            lin <- c(lin, getLineage(ped, founder, parents[1]))
        }
        if (parents[2] %in% allDes)
        {
            lin <- c(lin, getLineage(ped, founder, parents[2]))
        }
        return (lin)
    }
    else if (founder %in% parents)
    {
        return (parents)
    }
    else
    {
        return (c())
    }
}

getDistance <- function(ped, sub1, sub2)
{

}

getCommonAncestor <- function(ped, sub1, sub2)
{
    founders <- getFounders(ped)
    topFounders <- getTopFounders(ped)
    
    valid <- sapply(topFounders, function(x)
        {
            return (sub1 %in% getDescendants(ped, x)[1,]
                && sub2 %in% getDescendants(ped, x)[1,])
        })

    topFounders <- topFounders[valid]
    if (length(topFounders) == 0) return (NA)
    n <- sapply(topFounders, function(x) length(getDescendants(ped, x)[1,]))
    topFounder <- topFounders[n==max(n)][1]

    path1 <- getIntermediateAncestors(ped, topFounder, sub1)
    path2 <- getIntermediateAncestors(ped, topFounder, sub2)

    overlap <- c(path1, path2)
    overlap <- overlap[duplicated(overlap)]
    overlap <- overlap[!overlap %in% founders]

    if (length(overlap) == 0) return (topFounder)
    else return (overlap[1])
}

getBranchingPoints <- function(ped)
{
    finalDescendants <- combn(getFinalDescendants(ped), 2)
    branchingPoints <- apply(finalDescendants, 2, function(x)
        getCommonAncestor(ped, x[1], x[2]))

    branchingPoints <- unique(branchingPoints)
    founders <- getFounders(ped)
    return (branchingPoints[!(branchingPoints %in% founders)])
}

validPedigree <- function(ped)
{
    return (TRUE)
}

pedToDAG <- function(ped)
{
    # list used to generator DAG
    edgeList <- list()
    index <- 1

    # process nodes in order
    for (i in 1:getSize(ped))
    {
        edgeList[[index]] <- i
        index <- index + 1
    }

    # find edges in graph
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
