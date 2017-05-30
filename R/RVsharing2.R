# check if pedigree is valid for RVsharing
processPedigree <- function(ped)
{
    ped$id <- 1:length(ped$id)
    parents <- sapply(ped$id, function(i) c(ped$findex[i], ped$mindex[i]))
    founders <- which(parents[1,] == 0)
    affected <- which(ped$affected == 1)

    if (sum(affected %in% founders) > 0)
        stop('some founders are affected')
    if (sum(ped$affected==1) < 2)
        stop('need at least 2 affected subjects')

    return(list('ped'=ped, 'parents'=parents, 'founders'=founders,
        'affected'=affected))
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
    net <- compile(net)#, root=as.character(which(ped$affected==1)))
    return(propagate(net))
}

# prob of event that all marginal nodes are 0, and event all are 1
marginalProb <- function(net, marginalNodes)
{
    p0 <- p1 <- 1
    net0 <- net1 <- net
    for (n in as.character(marginalNodes))
    {
        if (p0 > 0) # prevents conditioning on zero prob events
        {
            p0 <- p0 * unname(querygrain(net0, n)[[1]][1])
            net0 <- setEvidence(net0, n, '0')
        }
        if (p1 > 0)
        {
            p1 <- p1 * unname(querygrain(net1, n)[[1]][2])
            net1 <- setEvidence(net1, n, '1')
        }
    }
    return(c(p0,p1))
}

oneFounderSharingProb <- function(net, founders, affected)
{
    # each fraction component of probability
    numer <- denom <- 0

    # set all founders to 0 (no variant)
    net <- setEvidence(net, as.character(founders),
        rep('0', length(founders)))

    # sum over probs, conditioning on each founder introducing variant
    for (f in founders)
    {
        # condition on founder and calculate distribution
        net <- retractEvidence(net, as.character(f))
        net <- setEvidence(net, as.character(f), '1')
        prob <- marginalProb(net, affected)

        # sum relevant probability       
        numer <- numer + prob[2]
        denom <- denom + 1 - prob[1]

        # reset founder
        net <- retractEvidence(net, as.character(f))
        net <- setEvidence(net, as.character(f), '0')
    }
    return(numer/denom)
}

pFU <- function(nFounders, theta, order=2)
{
    a <- (2*nFounders):(2*nFounders - order)
    dist <- c(1, theta, theta^2/2, theta^3/6, theta^4/24, theta^5/120)
    return(weighted.mean(2/nFounders - 2/a, dist[1:(order+1)]))
}

twoFounderSharingProb <- function(net, kinshipCoeff, founders,
affected, nSimulations)
{
    numer1 <- numer2 <- denom1 <- denom2 <- 0
    processed <- c()
    for (f1 in founders)
    {
        for (f2 in setdiff(founders, processed))
        {

        }
        processed <- c(processed, f1)
    }
}

#' @export
RVsharing2 <- function(ped, alleleFreq, kinshipCoeff, nSimulations)
{
    # pre-process pedigree
    ped <- processPedigree(ped)

    if (!missing(nSimulations))
        return(monteCarloSharingProb(ped, alleleFreq, kinshipCoeff, nSimulations))

    # calculate prior distribution based on allele frequency
    ifelse(missing(alleleFreq), p <- 0.5, p <- alleleFreq)
    prior <- c((1-p)^2, 2*p*(1-p), p^2)

    # create bayesian network from the pedigree
    net <- createNetwork(ped$ped, ped$parents, ped$founders, prior)

    # calculate sharing prob with appropiate method
    if (!missing(alleleFreq))
    {
        prob <- marginalProb(net, ped$affected)
        return(prob[2] / (1 - prob[1]))
    }        
    else if (!missing(kinshipCoeff))
    {
        # ... twoFounderSharingProb
    }
    else
    {
        return(oneFounderSharingProb(net, ped$founders, ped$affected))
    }
}



