context('PMF operations')

test_that('root pmf\'s correctly built',
{
    data(sampleGraphs)

    fixedVars <- matrix(0, nrow=2,ncol=0)
    prior <- lapply(1:10, function(x) c(0.25, 0.5, 0.25))
    pmf <- createPmf(sampleDAG, 1, fixedVars, prior)

    expect_equal(pmf[['vars']], c(1))
    expect_equal(pmf[['prob']], c(0.25, 0.5, 0.25))
})

test_that('non-root pmf\'s correctly built',
{
    data(sampleGraphs)

    fixedVars <- matrix(0, nrow=2,ncol=0)
    pmf <- createPmf(sampleDAG, 11, fixedVars)

    expect_equal(pmf[['vars']], c(11,8,13))
    expect_equal(dim(pmf[['prob']]), c(3,3,3))
    expect_equal(pmf[['prob']][1,1,1], 1)
    expect_equal(pmf[['prob']][1,3,3], 0)
    expect_equal(pmf[['prob']][1,2,2], 0.25)
    expect_equal(pmf[['prob']][2,2,2], 0.5)
    expect_equal(pmf[['prob']][3,2,2], 0.25)
})

test_that('non-root pmf\'s correctly built - account for fixed variables',
{
    data(sampleGraphs)

    fixedVars <- matrix(0, nrow=2,ncol=2)
    fixedVars[1,] <- 1:2
    pmf4 <- createPmf(sampleDAG, 4, fixedVars)
    pmf7 <- createPmf(sampleDAG, 7, fixedVars)
    pmf8 <- createPmf(sampleDAG, 8, fixedVars)

    expect_equal(pmf4[['vars']], c(4))
    expect_equal(pmf4[['prob']], c(1,0,0))

    expect_equal(pmf7[['vars']], c(7,3))
    expect_equal(dim(pmf7[['prob']]), c(3,3))
    expect_equal(pmf7[['prob']], matrix(c(1,0,0,0.5,0.5,0,0,1,0), ncol=3))

    expect_equal(pmf8[['vars']], c(8,4,5))
    expect_equal(dim(pmf8[['prob']]), c(3,3,3))
    expect_equal(pmf8[['prob']][1,1,1], 1)
    expect_equal(pmf8[['prob']][1,3,3], 0)
    expect_equal(pmf8[['prob']][1,2,2], 0.25)
    expect_equal(pmf8[['prob']][2,2,2], 0.5)
    expect_equal(pmf8[['prob']][3,2,2], 0.25)
})

test_that('probability accurately retrived from pmf',
{
    data(sampleGraphs)

    fixedVars <- matrix(0, nrow=2,ncol=2)
    fixedVars[1,] <- 1:2
    pmf4 <- createPmf(sampleDAG, 4, fixedVars)
    pmf7 <- createPmf(sampleDAG, 7, fixedVars)
    pmf8 <- createPmf(sampleDAG, 8, fixedVars)

    expect_equal(getProb(pmf4, c(4,1,2), c(1,0,0)),   0)
    expect_equal(getProb(pmf7, c(3,7),   c(1,0)),   0.5)
    expect_equal(getProb(pmf8, c(4,5,8), c(1,1,1)), 0.5)

    expect_error(getProb(pmf8, c(4,5), c(1,0)))
})

test_that('pmf\'s are multiplied correctly',
{
    data(sampleGraphs)

    fixedVars <- matrix(0, nrow=2,ncol=2)
    fixedVars[1,] <- 1:2
    pmf4 <- createPmf(sampleDAG, 4, fixedVars)
    pmf7 <- createPmf(sampleDAG, 7, fixedVars)
    pmf8 <- createPmf(sampleDAG, 8, fixedVars)
        
    pmf47 <- multiplyPmf(pmf4, pmf7)
    pmf48 <- multiplyPmf(pmf4, pmf8)
    pmf78 <- multiplyPmf(pmf7, pmf8)

    expect_equal(pmf47[['vars']], c(4,7,3))
    expect_equal(dim(pmf47[['prob']]), c(3,3,3))

    expect_equal(pmf48[['vars']], c(4,8,5))
    expect_equal(dim(pmf48[['prob']]), c(3,3,3))

    expect_equal(pmf78[['vars']], c(7,3,8,4,5))
    expect_equal(dim(pmf78[['prob']]), c(3,3,3,3,3))

    vars <- c(3,4,5,7,8)
    args <- expand.grid(lapply(1:5, function(x) 0:2))
    args <- as.matrix(unname(args))
    for (a in 1:nrow(args))
    {
        arg <- args[a,]
        expect_equal(getProb(pmf78, vars, arg),
        getProb(pmf7, vars, arg) * getProb(pmf8, vars, arg))
    }
})

test_that('variables are summed out correctly',
{
    data(sampleGraphs)

    fixedVars <- matrix(0, nrow=2,ncol=2)
    fixedVars[1,] <- 1:2
    pmf7 <- createPmf(sampleDAG, 7, fixedVars)
    pmf8 <- createPmf(sampleDAG, 8, fixedVars)
    pmf78 <- multiplyPmf(pmf7, pmf8)

    pmf78 <- sumOutVariable(pmf78, 3)
    pmf78 <- sumOutVariable(pmf78, 4)
    pmf78 <- sumOutVariable(pmf78, 5)
    pmf78 <- sumOutVariable(pmf78, 7)

    expect_equal(pmf78[['vars']], c(8))
    expect_equal(pmf78[['prob']], c(6.75,13.50,6.75))
})
