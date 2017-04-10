context('Pedigree Handling')

test_that('test basic pedigree operations',
{
    data(testPed)
    fam <- testPedList[[5]]

    expect_equal(getSize(fam), 27)

    expect_equal(getOffspring(fam, 19), 21)
    expect_equal(getOffspring(fam, 26), 27)
    expect_equal(getOffspring(fam, 27), integer(0))

    expect_equal(getParents(fam, 20), c(15, 16))
    expect_equal(getParents(fam, 11), c(7, 8))
    expect_equal(getParents(fam, 1), c(0, 0))

    expect_equal(getSpouses(fam, 3), 4)
    expect_equal(getSpouses(fam, 24), 23)
    expect_equal(getSpouses(fam, 27), 0)

    expect_equal(getInLaws(fam, 1), c(0, 0))
    expect_equal(getInLaws(fam, 16), c(11, 12))
    expect_equal(getInLaws(fam, 26), c(23, 24))
    
    expect_equal(getFounders(fam), c(1,2,4,6,8,10,12,14,16,19,22,24,26))
    expect_equal(getFinalDescendants(fam), c(17,20,27))

    expect_equal(getNextGeneration(fam, 1), c(3, 5))
    expect_equal(getNextGeneration(fam, c(3,5)), c(7, 9))
    expect_equal(getNextGeneration(fam, c(11,12)), c(15, 18))
    
    expect_equal(getDescendants(fam, 19)[1,], c(21, 23, 25, 27))
    expect_equal(getDescendants(fam, 19)[2,], c(1, 2, 3, 4))
})


