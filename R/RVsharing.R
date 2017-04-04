calculateSharingProb <- function(ped)
{
    affected <- which(ped$affected == 1) 

    founders <- getFounders(ped)   

    numer <- 0
    denom <- 0
    
    for (f in founders)
    {
        # get all descendants
        descendants <- getDescendants(ped, f)
        if (length(descendants) == 0)
        {
            stop('founders have no decendants')
        }

        # get distances of sequenced/affected descendants
        affectedCols <- which(descendants[1,] %in% affected)
        distances <- descendants[2, affectedCols]

        # calculate denominator
        denom <- denom + 1 - prod(1 - 0.5 ^ distances)

        # only add to numerator if founder is common
        if (all(affected %in% descendants[1,]))
        {
            numer <- numer + 0.5 ^ sum(distances)
        }
    }
    return (numer / denom)
}
