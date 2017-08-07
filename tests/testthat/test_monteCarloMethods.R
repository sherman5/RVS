context('monte carlo methods')

test_that('pedigree simulation',
{
    data(samplePedigrees)
    ped <- samplePedigrees$firstCousinPair
    procPed <- processPedigree(ped)

    states <- c(2,2,NA,2,NA,2,NA,NA)
    states <- simulatePedigree(procPed, states)
    expect_true(all(states==2))

    states <- c(0,0,NA,0,NA,0,NA,NA)
    states <- simulatePedigree(procPed, states)
    expect_true(all(states==0))

    states <- c(2,2,NA,0,NA,0,NA,NA)
    states <- simulatePedigree(procPed, states)
    expect_equal(states[procPed$carriers], c(1,1))
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
