getNonRootPmf <- function(var, parents, fixedVars)
{
    # create pmf structure
    pmf <- list()
    pmf[['vars']] <- c(var, parents)
    pmf[['prob']] <- array(0, c(3,3,3))

    # populate pmf probability
    args <- expand.grid(p1=0:2, p2=0:2)
    pmf[['prob']][1,,] <- with(args, (2-p1)*(2-p2)/4)
    pmf[['prob']][2,,] <- with(args, (p1+p2-p1*p2)/2)
    pmf[['prob']][3,,] <- with(args, p1*p2/4)

    # plug in fixed variable values
    pmf[['vars']] <- setdiff(pmf[['vars']], fixedVars[1,])
    i1 <- i2 <- i3 <- 1:3

    if (var %in% fixedVars[1,])
        i1 <- fixedVars[2,which(fixedVars[1,]==var)] + 1
    if (parents[1] %in% fixedVars[1,])
        i2 <- fixedVars[2,which(fixedVars[1,]==parents[1])] + 1
    if (parents[2] %in% fixedVars[1,])
        i3 <- fixedVars[2,which(fixedVars[1,]==parents[2])] + 1

    # prune pmf and return result
    pmf[['prob']] <- pmf[['prob']][i1,i2,i3]
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

multiplyPmf <- function(pmf1, pmf2)
{
    # create new pmf
    newPmf <- list()
    newPmf[['vars']] <- unique(c(pmf1[['vars']], pmf2[['vars']]))
    newPmf[['prob']] <- array(0, rep(3,length(newPmf[['vars']])))

    # find all variable combinations needed for new pmf
    args <- expand.grid(lapply(1:length(newPmf[['vars']]), function(x) 0:2))
    args <- as.matrix(unname(args))
    for (a in 1:nrow(args))
    {
        arg <- args[a,]
        newPmf[['prob']][matrix(arg+1, ncol=length(arg))] <-
            getProb(pmf1, newPmf[['vars']], arg) * 
            getProb(pmf2, newPmf[['vars']], arg)
    }
    return(newPmf)
}

createLocalPmf <- function(graph, vars, prior, fixedVars)
{
    # start with constant pmf
    pmf <- list()
    pmf[['prob']] <- c(1)

    # add each variable's indiviudal pmf
    for (v in vars)
    {
        parents <- as.numeric(parents(v,graph))
        if (length(parents) == 0) # founder
            indivPmf <- list('vars'=v, 'prob'= prior[[v]])
        else # non-founder
            indivPmf <- getNonRootPmf(v, parents, fixedVars)

        # update main pmf
        pmf <- multiplyPmf(pmf, indivPmf)
    }
    return(pmf)
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
