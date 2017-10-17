context('pedigree methods')

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
    expect_equal(procPed$carriers, c(7,8))
    expect_equal(procPed$id, 1:8)
    expect_equal(procPed$finalDescendants, c(7,8))
    
    ped$affected <- rep(0, procPed$size)
    expect_error(processPedigree(ped), 'need at least 2 affected subjects')
})

test_that('pedigree calculations',
{
    data(samplePedigrees)
    ped <- samplePedigrees$firstCousinPair
    procPed <- processPedigree(ped)

    # descendants
    expect_true(isDescendant(procPed, 1, 7))
    expect_true(isDescendant(procPed, 2, 8))
    expect_true(!isDescendant(procPed, 2, 4))
    expect_true(!isDescendant(procPed, 3, 5))
    expect_true(!isDescendant(procPed, 3, 8))
    
    # mating
    expect_true(areMating(procPed, 1, 2))
    expect_true(areMating(procPed, 3, 4))
    expect_true(!areMating(procPed, 2, 4))
    expect_true(!areMating(procPed, 3, 5))
    expect_true(!areMating(procPed, 6, 2))

    # distance
    expect_equal(ancestorDistance(procPed, 1, 3), 1)
    expect_equal(ancestorDistance(procPed, 2, 7), 2)
    expect_equal(ancestorDistance(procPed, 1, 5), 1)
    expect_equal(ancestorDistance(procPed, 2, 8), 2)
    expect_error(ancestorDistance(procPed, 1, 6))
    expect_error(ancestorDistance(procPed, 5, 6))
    expect_error(ancestorDistance(procPed, 8, 6))
    expect_error(ancestorDistance(procPed, 1, 2))
})

test_that('kinship coefficient estimation',
{
    data(samplePedigrees)
    
    coefMatrix <- ComputeKinshipPropCoef(samplePedigrees$firstCousinPair)
    expect_equal(unname(coefMatrix[1,]), c(   NA, 0.125))
    expect_equal(unname(coefMatrix[2,]), c(0.125,    NA))
})


