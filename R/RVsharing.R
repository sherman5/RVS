# construct bayesian network using gRain package
buildBayesNet <- function(ped, prior)
{
    # get list of all parents, determine founders
    ped$id <- 1:length(ped$id)
    parents <- sapply(ped$id, function(i) c(ped$findex[i], ped$mindex[i]))
    founders <- which(parents[1,] == 0)

    # process founders
    founderNodes <- lapply(founders, cptable, values=prior, levels=0:2)

    # conditional prob table based on mendel's laws
    tbl <- array(0, c(3,3,3))
    args <- expand.grid(p1=0:2, p2=0:2)
    tbl[1,,] <- with(args, (2-p1) * (2-p2) / 4)
    tbl[2,,] <- with(args, (p1 + p2 - p1*p2) / 2)
    tbl[3,,] <- with(args, p1 * p2 / 4)

    # process non-founders
    nonFounderNodes <- lapply(setdiff(ped$id, founders),
        function(nf) {cptable(c(nf, parents[1,nf], parents[2,nf]),
        values=c(tbl), levels=0:2)})

    # create bayesian network
    return(grain(compileCPT(c(founderNodes, nonFounderNodes))))
}

# calculate standard sharing probability on the assumption of 1 founder
# introducing the variant
sharingProb <- function(ped)
{
    # process pedigree
    ped$id <- 1:length(ped$id)
    affected <- which(ped$affected == 1)
    parents <- sapply(ped$id, function(i) c(ped$findex[i], ped$mindex[i]))
    founders <- which(parents[1,] == 0)

    # construct bayesian network, TODO: don't set root if # affected large
    bayesNet <- buildBayesNet(ped, c(1,0,0))
    bayesNet <- compile(bayesNet, root=as.character(affected))
    bayesNet <- propagate(bayesNet)

    # condition on each founder introducing the variant
    numer <- 0
    denom <- 0
    numAff <- length(affected)
    for (f in founders)
    {
        # condition on founder and calculate distribution
        bayesNet <- setEvidence(bayesNet, as.character(f), '1')
        query <- querygrain(bayesNet, as.character(affected), type='joint')

        # sum relevant probability       
        numer <- numer + query[matrix(rep(2,numAff), ncol=numAff)]
        denom <- denom + 1 - query[matrix(rep(1,numAff), ncol=numAff)]

        # reset founder
        bayesNet <- retractEvidence(bayesNet, as.character(f))
    }
    return(numer/denom)
}

createNetwork <- function(ped, alleleFreq)
{

}

evaluateNetwork <- function(net, nSimulations)
{


}

#' @export
RVsharing <- function(ped, alleleFreq, nSimulations)
{
    if (!missing(alleleFreq))
    {

    }
    else if (!missing(nSimulations))
    {

    }
    else
    {

    }
    return(sharingProb)
}



