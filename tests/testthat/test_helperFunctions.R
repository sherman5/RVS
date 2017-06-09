context('helper functions')

test_that('pre-process pedigree',
{
    data(samplePedigrees)

    ped <- samplePedigrees$firstCousinPair
    procPed <- processPedigree(ped)

    expect_equal(procPed$ped, ped)
    expect_equal(procPed$parents[,1], c(0,0))
    expect_equal(procPed$parents[,2], c(0,0))
    expect_equal(procPed$parents[,3], c(1,2))
    expect_equal(procPed$parents[,5], c(1,2))
    expect_equal(procPed$parents[,7], c(3,4))
    expect_equal(procPed$parents[,8], c(5,6))
    expect_equal(procPed$founders, c(1,2,4,6))
    expect_equal(procPed$affected, c(7,8))
    expect_equal(procPed$size, 8)
    
    ped$affected[1] <- 1
    expect_error(RVsharing2(ped), 'some founders are affected')

    ped$affected <- rep(0, procPed$size)
    expect_error(RVsharing2(ped), 'need at least 2 affected subjects')
})

test_that('pedigree simulation',
{
    data(samplePedigrees)
    ped <- samplePedigrees$firstCousinPair
    procPed <- processPedigree(ped)

    states <- c(2,2,NA,2,NA,2,NA,NA)
    simPed <- simulatePedigree(procPed, states)
    expect_equal(simPed, c(2,2))

    states <- c(0,0,NA,0,NA,0,NA,NA)
    simPed <- simulatePedigree(procPed, states)
    expect_equal(simPed, c(0,0))

    states <- c(2,2,NA,0,NA,0,NA,NA)
    simPed <- simulatePedigree(procPed, states)
    expect_equal(simPed, c(1,1))
})

test_that('full monte carlo simulation',
{
    data(samplePedigrees)
    ped <- samplePedigrees$firstCousinPair
    procPed <- processPedigree(ped)

    founderFunc <- function(n) rep(0,n)
    simProb <- runMonteCarlo(procPed, founderFunc, 10)
    expect_equal(simProb, NaN)

    founderFunc <- function(n) rep(2,n)
    simProb <- runMonteCarlo(procPed, founderFunc, 10)
    expect_equal(simProb, 1)
})

test_that('building a bayesian network',
{
    data(samplePedigrees)
    ped <- samplePedigrees$firstCousinPair
    procPed <- processPedigree(ped)

    net <- createNetwork(procPed)
    expect_equal(length(net$universe$nodes), procPed$size)

    net <- createNetwork(procPed, c(1,0,0))
    expect_equal(unname(gRain::querygrain(net, '7')[[1]][1]), 1)
    expect_equal(unname(gRain::querygrain(net, '8')[[1]][1]), 1)
    expect_equal(unname(gRain::querygrain(net, '7')[[1]][3]), 0)
    expect_equal(unname(gRain::querygrain(net, '8')[[1]][3]), 0)
    
    net <- createNetwork(procPed, c(0,0,1))
    expect_equal(unname(gRain::querygrain(net, '7')[[1]][3]), 1)
    expect_equal(unname(gRain::querygrain(net, '8')[[1]][3]), 1)
})

test_that('marginal probability calculation',
{
    data(samplePedigrees)
    ped <- samplePedigrees$firstCousinPair
    procPed <- processPedigree(ped)
    net <- createNetwork(procPed)

    expect_equal(marginalProb(net, list('7'=0)), 1/4)
    expect_equal(marginalProb(net, list('8'=c(1,2))), 3/4)
    expect_equal(marginalProb(net, list('4'=0,'6'=0)), 1/16)
    expect_true(all(marginalProb(net, c(7,8)) >= marginalProb(net, c(4,6))))
})
