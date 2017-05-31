context('helper functions')

test_that('pre-process pedigree',
{
    data(samplePedigrees)

    ped <- samplePedigrees$firstCousinPair
    procPed <- processPedigree(ped)
    # should run fine, inspect output
    
    ped$affected[1] <- 1
    # expect error where founder is affected

    ped$affected <- rep(0, procPed$size)
    # expect error where not enough subjects are affected
})

test_that('pedigree simulation',
{
    data(samplePedigrees)

    
})

test_that('full monte carlo simulation',
{

})

test_that('building a bayesian network',
{

})

test_that('marginal probability calculation',
{

})
