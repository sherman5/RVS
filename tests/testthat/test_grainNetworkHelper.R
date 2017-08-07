context('grain network helper functions')

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

    expect_equal(marginalProb(net, list('1'=0)), 1/4)
    expect_equal(marginalProb(net, list('1'=1)), 1/2)
    expect_equal(marginalProb(net, list('1'=2)), 1/4)

    expect_equal(marginalProb(net, list('3'=0)), 1/4)
    expect_equal(marginalProb(net, list('3'=1)), 1/2)
    expect_equal(marginalProb(net, list('3'=2)), 1/4)

    expect_equal(marginalProb(net, list('7'=0)), 1/4)
    expect_equal(marginalProb(net, list('7'=1)), 1/2)
    expect_equal(marginalProb(net, list('7'=2)), 1/4)

    expect_equal(marginalProb(net, list('6'=0, '8'=c(1,2))), 1/8)
    expect_equal(marginalProb(net, list('6'=2, '8'=c(1,2))), 1/4)
    expect_equal(marginalProb(net, list('6'=c(0,1), '4'=c(1,2))), 9/16)
})
