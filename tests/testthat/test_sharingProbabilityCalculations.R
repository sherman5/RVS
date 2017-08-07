context('standard sharing probability calculations')

test_that('exact sharing probability',
{
    data(samplePedigrees)
       
    expect_equal(RVsharing(samplePedigrees$firstCousinPair,
        alleleFreq=0.01), 0.076, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$secondCousinPair,
        alleleFreq=0.01), 0.025, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$thirdCousinPair,
        alleleFreq=0.01), 0.014, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$firstCousinTriple,
        alleleFreq=0.01), 0.014, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$secondCousinTriple,
        alleleFreq=0.01), 0.002, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$secondCousinPairWithLoop,
        alleleFreq=0.01), 0.058, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$firstCousinInbreeding,
        alleleFreq=0.01), 0.122, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$twoGenerationsInbreeding,
        alleleFreq=0.01), 0.029, tol=1e-3)
})

test_that('one founder sharing probability',
{
    data(samplePedigrees)

    expect_equal(RVsharing(samplePedigrees$firstCousinPair),
        1/15)
    expect_equal(RVsharing(samplePedigrees$secondCousinPair),
        1/63)
    expect_equal(RVsharing(samplePedigrees$thirdCousinPair),
        1/255)
    expect_equal(RVsharing(samplePedigrees$firstCousinTriple),
        1/85)
    expect_equal(RVsharing(samplePedigrees$secondCousinTriple),
        1/745)
    expect_equal(RVsharing(samplePedigrees$secondCousinPairWithLoop),
        6/122)
    expect_equal(RVsharing(samplePedigrees$firstCousinInbreeding),
        201/1783)
    expect_equal(RVsharing(samplePedigrees$twoGenerationsInbreeding),
        37/1975)
})


