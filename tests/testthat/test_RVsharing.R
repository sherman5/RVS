context('Top level function')

test_that('standard sharing probabilities',
{
    data(samplePedigrees)

    expect_equal(RVsharing(samplePedigrees$firstCousinPair), 1/15)
    expect_equal(RVsharing(samplePedigrees$secondCousinPair), 1/63)
    expect_equal(RVsharing(samplePedigrees$thirdCousinPair), 1/255)
    expect_equal(RVsharing(samplePedigrees$firstCousinTriple), 1/85)
    expect_equal(RVsharing(samplePedigrees$secondCousinTriple), 1/745)
    expect_equal(RVsharing(samplePedigrees$secondCousinPairWithLoop), 6/122)
    expect_equal(RVsharing(samplePedigrees$firstCousinInbreeding), 180/1783)
    expect_equal(RVsharing(samplePedigrees$twoGenerationsInbreeding), 34/1975)
})



