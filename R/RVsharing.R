# construct bayesian network using gRain package
buildBayesNet <- function(ped)
{
    # get list of all parents, determine founders
    ped$id <- 1:length(ped$id)
    parents <- sapply(ped$id, function(i) c(ped$findex[i], ped$mindex[i]))
    founders <- which(parents[1,] == 0)

    # process founders
    founderNodes <- lapply(founders, cptable, values=c(1,2,1),
        levels=0:2)

    # conditional prob table based on mendel's laws
    tbl <- array(0, c(3,3,3))
    args <- expand.grid(p1=0:2, p2=0:2)
    tbl[1,,] <- with(args, (2-p1)*(2-p2)/4)
    tbl[2,,] <- with(args, (p1+p2-p1*p2)/2)
    tbl[3,,] <- with(args, p1*p2/4)

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

    # construct bayesian network
    bayesNet <- buildBayesNet(ped)
    bayesNet <- compile(bayesNet, root=as.character(affected))
    bayesNet <- propagate(bayesNet)
    bayesNet <- setEvidence(bayesNet, as.character(founders),
        rep('0', length(founders)))

    # condition on each founder introducing the variant
    numer <- 0
    denom <- 0
    for (f in founders)
    {
        # condition on founder and calculate distribution
        bayesNet <- retractEvidence(bayesNet, as.character(f))
        bayesNet <- setEvidence(bayesNet, as.character(f), '1')
        query <- querygrain(bayesNet, as.character(affected), type='joint')

        # sum relevant probability       
        n <- length(affected)
        numer <- numer + query[matrix(rep(2, n), ncol=n)]
        denom <- denom + 1 - query[matrix(rep(1, n), ncol=n)]

        # reset founder
        bayesNet <- retractEvidence(bayesNet, as.character(f))
        bayesNet <- setEvidence(bayesNet, as.character(f), '0')
    }
    return(numer/denom)
}


