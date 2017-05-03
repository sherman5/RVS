context('PMF operations')

test_that('non-root pmf\'s correctly built',
{
    fixedVars <- matrix(0, nrow=2,ncol=1)
    pmf <- getNonRootPmf(1, c(2,3), fixedVars)

    expect_equal(dim(pmf[['prob']]), c(3,3,3))
    expect_equal(pmf[['prob']][1,1,1], 1)
    expect_equal(pmf[['prob']][1,3,3], 0)
    expect_equal(pmf[['prob']][1,2,2], 0.25)
    expect_equal(pmf[['prob']][2,2,2], 0.5)
    expect_equal(pmf[['prob']][3,2,2], 0.25)
})

test_that('non-root pmf\'s correctly built - account for fixed variables',
{
    fixedVars <- matrix(0,2,1)
    fixedVars[1,] <- c(3)
    pmf <- getNonRootPmf(1, c(2,3), fixedVars)
    expect_equal(dim(pmf[['prob']]), c(3,3))
    expect_equal(pmf[['prob']], matrix(c(1,0,0,0.5,0.5,0,0,1,0),3,3))

    fixedVars <- matrix(0,2,2)
    fixedVars[1,] <- c(2,3)
    pmf <- getNonRootPmf(1, c(2,3), fixedVars)
    expect_equal(pmf[['prob']], c(1,0,0))

    fixedVars <- matrix(0,2,3)
    fixedVars[1,] <- c(1,2,3)
    pmf <- getNonRootPmf(1, c(2,3), fixedVars)
    expect_equal(pmf[['prob']], 1)
})

test_that('probability accurately retrived from pmf',
{
    fixedVars <- matrix(0, nrow=2,ncol=1)
    pmf <- getNonRootPmf(1, c(2,3), fixedVars)

    vars <- c(1,2,3)
    expect_equal(getProb(pmf, vars, c(0,0,0)), 1)
    expect_equal(getProb(pmf, vars, c(0,2,2)), 0)
    expect_equal(getProb(pmf, vars, c(2,2,2)), 1)
    expect_equal(getProb(pmf, vars, c(0,1,1)), 0.25)
    expect_equal(getProb(pmf, vars, c(1,1,1)), 0.5)
    expect_equal(getProb(pmf, vars, c(2,1,1)), 0.25)
})

test_that('pmf\'s are multiplied correctly',
{
    fixedVars <- matrix(0, nrow=2,ncol=1)
    pmf1 <- getNonRootPmf(1, c(2,3), fixedVars)
    pmf2 <- getNonRootPmf(3, c(4,5), fixedVars)
        
    prodPmf <- multiplyPmf(pmf1, pmf2)
    expect_equal(prodPmf[['vars']], 1:5)
    expect_equal(dim(prodPmf[['prob']]), c(3,3,3,3,3))

    args <- expand.grid(lapply(1:5, function(x) 0:2))
    args <- as.matrix(unname(args))
    for (a in 1:nrow(args))
    {
        arg <- args[a,]
        expect_equal(getProb(prodPmf, 1:5, arg),
            getProb(pmf1, 1:5, arg) * getProb(pmf2, 1:5, arg))
    }
})

test_that('local pmf\'s correctly built',
{
    graph <- dag(c(1,2,3), c(3,4,5), c(4,6,7), c(10,6,7), c(9,8,10))
    fixedVars <- matrix(0,2,2)
    fixedVars[1,] <- c(1,9)
    fixedVars[2,] <- c(0,0)
    prior <- lapply(1:10, function(x) c(0.25,0.5,0.25))
    pmf <- createLocalPmf(graph, c(1,2,3), prior, fixedVars)

    f1 <- function(x,a,b) c((2-a)*(2-b)/4,(a+b-a*b)/2,a*b/4)[x+1]
    f2 <- function(x) c(0.25,0.5,0.25)[x+1]
    args <- expand.grid(lapply(1:4, function(x) 0:2))
    args <- as.matrix(unname(args))
    for (a in 1:nrow(args))
    {
        arg <- args[a,]
        expect_equal(getProb(pmf, 2:5, arg), f1(0,arg[1],arg[2]) * 
            f1(arg[2],arg[3],arg[4]) * f2(arg[1]))
    }
})

test_that('variables are summed out correctly',
{
    graph <- dag(c(1,2,3), c(3,4,5), c(4,6,7), c(10,6,7), c(9,8,10))
    fixedVars <- matrix(0,2,2)
    fixedVars[1,] <- c(1,9)
    fixedVars[2,] <- c(0,0)
    prior <- lapply(1:10, function(x) c(0.25,0.5,0.25))
    pmf <- createLocalPmf(graph, c(1,2,3), prior, fixedVars)

    pmf <- sumOutVariable(pmf,4)
    pmf <- sumOutVariable(pmf,5)

    expect_equal(pmf[['vars']], c(2,3))
    expect_equal(pmf[['prob']],
        matrix(c(0.5625,0.5625,0,0.5625,0.5625,0,0,0,0)),3,3)
})
