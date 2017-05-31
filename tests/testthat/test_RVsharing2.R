context('top level function')

test_that('standard sharing probabilities',
{
    data(samplePedigrees)

    expect_equal(RVsharing2(samplePedigrees$firstCousinPair), 1/15)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPair), 1/63)
    expect_equal(RVsharing2(samplePedigrees$thirdCousinPair), 1/255)
    expect_equal(RVsharing2(samplePedigrees$firstCousinTriple), 1/85)
    expect_equal(RVsharing2(samplePedigrees$secondCousinTriple), 1/745)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPairWithLoop),
        6/122)
    expect_equal(RVsharing2(samplePedigrees$firstCousinInbreeding),
        180/1783)
    expect_equal(RVsharing2(samplePedigrees$twoGenerationsInbreeding),
        34/1975)
})

test_that('monte carlo sharing probabilities',
{
    data(samplePedigrees)
       
    expect_equal(RVsharing2(samplePedigrees$firstCousinPair,
        nSimulations=1e4), 1/15, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPair,
        nSimulations=1e4), 1/63, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$thirdCousinPair,
        nSimulations=1e4), 1/255, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$firstCousinTriple,
        nSimulations=1e4), 1/85, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$secondCousinTriple,
        nSimulations=1e4), 1/745, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPairWithLoop,
        nSimulations=1e4), 6/122, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$firstCousinInbreeding,
        nSimulations=1e4), 180/1783, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$twoGenerationsInbreeding,
        nSimulations=1e4), 34/1975, tolerance=1e-2)

    expect_equal(RVsharing2(samplePedigrees$firstCousinPair,
        alleleFreq=0.05, nSimulations=1e4), 1/15, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPair,
        alleleFreq=0.05, nSimulations=1e4), 1/15, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$thirdCousinPair,
        alleleFreq=0.05, nSimulations=1e4), 1/15, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$firstCousinTriple,
        alleleFreq=0.05, nSimulations=1e4), 1/15, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$secondCousinTriple,
        alleleFreq=0.05, nSimulations=1e4), 1/15, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPairWithLoop,
        alleleFreq=0.05, nSimulations=1e4), 1/15, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$firstCousinInbreeding,
        alleleFreq=0.05, nSimulations=1e4), 1/15, tolerance=1e-2)
    expect_equal(RVsharing2(samplePedigrees$twoGenerationsInbreeding,
        alleleFreq=0.05, nSimulations=1e4), 1/15, tolerance=1e-2)
})

test_that('exact sharing probabilities',
{
    data(samplePedigrees)
       
    expect_equal(RVsharing2(samplePedigrees$firstCousinPair,
        alleleFreq=0.05), 1/15)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPair,
        alleleFreq=0.05), 1/15)
    expect_equal(RVsharing2(samplePedigrees$thirdCousinPair,
        alleleFreq=0.05), 1/15)
    expect_equal(RVsharing2(samplePedigrees$firstCousinTriple,
        alleleFreq=0.05), 1/15)
    expect_equal(RVsharing2(samplePedigrees$secondCousinTriple,
        alleleFreq=0.05), 1/15)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPairWithLoop,
        alleleFreq=0.05), 1/15)
    expect_equal(RVsharing2(samplePedigrees$firstCousinInbreeding,
        alleleFreq=0.05), 1/15)
    expect_equal(RVsharing2(samplePedigrees$twoGenerationsInbreeding,
        alleleFreq=0.05), 1/15)
})

