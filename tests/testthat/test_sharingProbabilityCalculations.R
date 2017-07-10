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

test_that('two founder sharing probability',
{
    data(samplePedigrees)

    phi <- 0.01

    expect_equal(RVsharing(samplePedigrees$firstCousinPair,
        kin=phi), 0.074, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$secondCousinPair,
        kin=phi), 0.0228, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$thirdCousinPair,
        kin=phi), 0.0106, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$firstCousinTriple,
        kin=phi), 0.0134, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$secondCousinTriple,
        kin=phi), 0.00174, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$secondCousinPairWithLoop,
        kin=phi), 0.056, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$firstCousinInbreeding,
        kin=phi), 0.119, tol=1e-3)
    expect_equal(RVsharing(samplePedigrees$twoGenerationsInbreeding,
        kin=phi), 0.0251, tol=1e-3)
})
