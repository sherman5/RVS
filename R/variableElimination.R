RVsharing <- function(ped)
{
    ped[['id']] <- 1:length(fam[['id']])
    affected <- which(ped$affected == 1)
    margNodes <- matrix(nrow=2, ncol=length(affected))
    margNodes[1,] <- affected
    prior <- lapply(1:getSize(fam), function(x) c(1,0,0))
    
    founders <- getFounders(ped)
    numer <- 0
    denom <- 0
    for (f in founders)
    {
        prior[[f]] <- c(0,1,0)
        margNodes[2,] <- rep(0, length(affected))
        margProb_0 <- marginalProbability(ped, margNodes, prior)
        margNodes[2,] <- rep(1, length(affected))
        margProb_1 <- marginalProbability(ped, margNodes, prior)  
        prior[[f]] <- c(1,0,0)

        numer <- numer + margProb_1
        denom <- denom + 1 - margProb_0
    }
    return(numer/denom)
}

marginalProbability <- function(ped, marginalVars, prior)
{
    # get graph order
    graph <- pedToDAG(ped)
    order <- as.numeric(topoSort(graph))
    order <- order[!order %in% marginalNodes[1,]]

    # partial sum pmf
    pmf <- list()
    pmf[['vars']] <- c()
    pmf[['prob']] <- array()

    # already processed nodes
    processed <- c()

    # sum over variables in order
    for (node in order)
    {
        # find all un-processed children for this node
        candidates <- c(node, getOffspring(ped, node))
        vars <- candidates[sapply(candidates, function(x)
            {sum(getParents(ped, x) %in% processed) == 0})]

        # find the parents that vars depend on, create local pmf
        varsDep <- unique(c(sapply(vars, getParents, ped=ped)))
        varsDep <- unlist(varsDep[localVarsDep != 0])
        localPmf <- createLocalPmf(ped, vars, varsDep, prior, marginalVars)

        # update old pmf with new local pmf
        pmf <- updatePmf(pmf, localPmf)

        # sum out node, update pmf, mark node as processed
        pmf <- sumOutVariable(pmf, node)
        processed <- c(processed, node)
    }
    return (pmf)
}

updatePmf <- function(pmf, newPmf)
{
    # create new pmf
    newPmf <- list()
    newPmf[['vars']] <- unique(c(pmf[['vars']], newPmf[['vars']]))
    newPmf[['prob']] <- array(0, rep(3,length(newVars)))

    # find all variable combinations needed for new pmf
    args <- expand.grid(lapply(1:length(newPmf[['vars']]), function(x) 0:2))
    args <- as.matrix(unname(args))
    for (a in 1:nrow(args))
    {
        arg <- args[a,]
        newPmf[['prob']][matrix(arg, ncol=length(arg))] <-
            getProb(pmf, newPmf[['vars']], arg) * 
            getProb(newPmf, newPmf[['vars']], arg)
    }
    return(newPmf)
}

sumOutVariable <- function(pmf, var)
{
    if (length(pmf[['vars']]) == 1)
    {
        pmf[['vars']] <- c()
        pmf[['prob']] <- sum(pmf[['prob']])
    }
    else
    {
        index <- which(pmf[['vars']] == var)
        perm <- 1:length(dim(pmf[['prob']]))
        perm <- perm[perm != index]
        perm <- c(index, perm)
        pmf[['prob']] <- colSums(aperm(pmf[['prob']], perm))
        pmf[['vars']] <- pmf[['vars']][pmf[['vars']] != var]
    }
    return(pmf)
}

getProb <- function(pmf, vars, vals)
{
    if (length(pmf[['vars']]) == 0)
    {
        return (pmf[['prob']])
    }
    else if (!all(pmf[['vars']] %in% vars))
    {
        stop('not enough variables supplied to pmf')
    }
    else
    {
        ind <- sapply(pmf[['vars']], function(x) which(vars==x))
        arrayIndex <- matrix(vals[ind] + 1, ncol=length(pmf[['vars']]))
        return (pmf[['prob']][arrayIndex])
    }
}

# only called on non-founders
getIndividualPmf <- function(var, parents, marginalVars)
{
    # create pmf structure
    pmf <- list()
    pmf[['vars']] <- c(var, parents)
    pmf[['prob']] <- array(0, c(3,3,3))

    # populate pmf args
    args <- expand.grid(x=0:2, a=0:2, b=0:2)
    pmf[['prob']][] <- with(args, c((2-a)*(2-b)/4, (a+b-a*b)/2, a*b/4)[x+1])   

    # plug in marginal variable values
    varIndex <- p1Index <- p2Index <- 1:3
    
    if (var %in% marginalVars[1,])
    {
        varIndex <- marginalVars[2,which(marginalVars[1,]==var)] + 1
        pmf[['vars']] <- pmf[['vars']][pmf[['vars']] != var]
    }
    if (parents[1] %in% marginalVars[1,])
    {
        p1Index <- marginalVars[2,which(marginalVars[1,]==parents[1])] + 1
        pmf[['vars']] <- pmf[['vars']][pmf[['vars']] != parents[1]]
    }
    if (parents[2] %in% marginalVars[1,])
    {
        p2Index <- marginalVars[2,which(marginalVars[1,]==parents[2])] + 1
        pmf[['vars']] <- pmf[['vars']][pmf[['vars']] != parents[2]]
    }

    pmf[['prob']] <- pmf[['prob']][varIndex, p1Index, p2Index]
    return(pmf)
}

createLocalPmf <- function(ped, vars, prior, marginalVars)
{
    pmf <- list()
    pmf[['prob']] <- c(1)

    for (v in vars)
    {
        # create PMF for this variable
        indivPmf <- list()
        if (sum(getParents(ped, v)) == 0)
        {
            indivPmf[['vars']] <- v
            indivPmf[['prob']] <- prior[[v]]
        }
        else
        {
            indivPmf <- getIndividualPmf(ped, v, prior, marginalVars)
        } 

        # update main pmf
        pmf <- updatePmf(pmf, indivPmf)
    }
    return(pmf)
}
