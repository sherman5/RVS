# check if pedigree is valid for RVsharing
validPedigree <- function(ped)
{
    if (sum(ped$affected==1) < 2) stop('need at least 2 affected subjects')
}

# create bayesian network from pedigree, use gRain package
createNetwork <- function(ped, parents, founders, prior)
{
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
    net <- grain(compileCPT(c(founderNodes, nonFounderNodes)))
    net <- compile(net, root=as.character(which(ped$affected==1)))
    return(propagate(net))
}

# calculate the joint-marginal distribution of specified nodes
marginalProb <- function(net, marginalNodes, nSimulations)
{
    # calculate exact distribution
    if (missing(nSimulations))
    {
        prob <- querygrain(net, as.character(marginalNodes), type='joint')
        p0 <-  prob[matrix(rep(1,length(marginalNodes)),
            ncol=length(marginalNodes))]
        p1 <-  prob[matrix(rep(2,length(marginalNodes)),
            ncol=length(marginalNodes))]
    }
    else # calculate distribution from monte carlo simulations
    {
        sim <- simulate(net, nsim=nSimulations)[,marginalNodes]
        sim <- matrix(as.numeric(as.matrix(sim)),ncol=length(marginalNodes))
        p0 <- sum(apply(sim, 1, function(r) all(r==0))) / nSimulations
        p1 <- sum(apply(sim, 1, function(r) all(r==1))) / nSimulations
    }
    return(c(p0,p1))
}

#' @export
RVsharing <- function(ped, alleleFreq, nSimulations)
{
    # check pedigree
    ped$id <- 1:length(ped$id)
    validPedigree(ped)
    parents <- sapply(ped$id, function(i) c(ped$findex[i], ped$mindex[i]))
    founders <- which(parents[1,] == 0)
    affected <- which(ped$affected == 1)
    nAff <- length(affected)
    
    # calculate prior distribution based on allele frequency
    ifelse(missing(alleleFreq), p <- 0.5, p <- alleleFreq)
    prior <- c((1-p)^2, 2*p*(1-p), p^2)

    # create bayesian network from the pedigree
    net <- createNetwork(ped, parents, founders, prior)

    # if allele frequency not given, assume 1 founder introduces variant
    if (missing(alleleFreq))
    {
        # each fraction component of probability
        numer <- 0
        denom <- 0

        # condition on each founder introducing the variant
        net <- setEvidence(net, as.character(founders),
            rep('0', length(founders)))
        for (f in founders)
        {
            # condition on founder and calculate distribution
            net <- retractEvidence(net, as.character(f))
            net <- setEvidence(net, as.character(f), '1')
            prob <- marginalProb(net, affected, nSimulations)

            # sum relevant probability       
            numer <- numer + prob[2]
            denom <- denom + 1 - prob[1]

            # reset founder
            net <- retractEvidence(net, as.character(f))
            net <- setEvidence(net, as.character(f), '0')
        }
        return(numer/denom)
    }
    else
    {
        prob <- marginalProb(net, affected, nSimulations)
        return(p[2] / (1 - p[1]))
    }
}



