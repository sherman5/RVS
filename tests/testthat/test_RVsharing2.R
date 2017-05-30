context('Top level function')

test_that('standard sharing probabilities',
{
    data(samplePedigrees)

    expect_equal(RVsharing2(samplePedigrees$firstCousinPair), 1/15)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPair), 1/63)
    expect_equal(RVsharing2(samplePedigrees$thirdCousinPair), 1/255)
    expect_equal(RVsharing2(samplePedigrees$firstCousinTriple), 1/85)
    expect_equal(RVsharing2(samplePedigrees$secondCousinTriple), 1/745)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPairWithLoop), 6/122)
    expect_equal(RVsharing2(samplePedigrees$firstCousinInbreeding), 180/1783)
    expect_equal(RVsharing2(samplePedigrees$twoGenerationsInbreeding), 34/1975)
})
