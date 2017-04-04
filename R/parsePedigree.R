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
    if (spouses[1] == 0) {return (c(0))}

    inLaws <- sapply(spouses, getParents, ped = ped)
    return (unique(c(inLaws)))
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

