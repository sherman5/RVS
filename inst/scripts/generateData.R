library(kinship2)

simPedigree <- function(numGenerations, numChildren, affRate)
{
    ped <- list()
    ped$id <- c(1)
    ped$findex <- c(0)
    ped$mindex <- c(0)
    ped$sex <- c(1)
    ped$affected <- c(0)
    class(ped) <- 'pedigree'

    moms <- c()
    curGen <- c(1)
    back <- function(v) v[length(v)]

    for (gen in 1:numGenerations)
    {
        prevGen <- curGen
        curGen <- c()
        for (f in prevGen)
        {
            children <- 1:numChildren() + back(ped$id)
            mom <- back(children) + 1

            curGen <- c(curGen, children)
            ped$id <- c(ped$id, children, mom)

            ped$findex <- c(ped$findex, rep(f, length(children)))
            ped$mindex <- c(ped$mindex, rep(mom, length(children)))
            ped$sex <- c(ped$sex, rep(1, length(children)))
            aff <- sample(0:1, length(children), replace=TRUE,
                prob = c(1-affRate, affRate))
            ped$affected <- c(ped$affected, aff)

            ped$findex <- c(ped$findex, 0)
            ped$mindex <- c(ped$mindex, 0)
            ped$sex <- c(ped$sex, 2)
            ped$affected <- c(ped$affected, 0)
        }
    }
    return(ped)
}

simInbredPedigree <- function(degree)
{
    ped <- list()
    ped$id <- c(1)
    ped$findex <- c(0)
    ped$mindex <- c(0)
    ped$sex <- c(1)
    ped$affected <- c(0)
    class(ped) <- 'pedigree'

}
