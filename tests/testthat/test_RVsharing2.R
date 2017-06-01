context('top level function')

test_that('standard sharing probabilities',
{
    data(samplePedigrees)

    expect_equal(RVsharing2(samplePedigrees$firstCousinPair),
        1/15)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPair),
        1/63)
    expect_equal(RVsharing2(samplePedigrees$thirdCousinPair),
        1/255)
    expect_equal(RVsharing2(samplePedigrees$firstCousinTriple),
        1/85)
    expect_equal(RVsharing2(samplePedigrees$secondCousinTriple),
        1/745)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPairWithLoop),
        6/122)
    expect_equal(RVsharing2(samplePedigrees$firstCousinInbreeding),
        201/1783)
    expect_equal(RVsharing2(samplePedigrees$twoGenerationsInbreeding),
        37/1975)
})

test_that('exact sharing probabilities',
{
    data(samplePedigrees)
       
    expect_equal(RVsharing2(samplePedigrees$firstCousinPair,
        alleleFreq=0.01), 0.076, tol=1e-3)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPair,
        alleleFreq=0.01), 0.025, tol=1e-3)
    expect_equal(RVsharing2(samplePedigrees$thirdCousinPair,
        alleleFreq=0.01), 0.014, tol=1e-3)
    expect_equal(RVsharing2(samplePedigrees$firstCousinTriple,
        alleleFreq=0.01), 0.014, tol=1e-3)
    expect_equal(RVsharing2(samplePedigrees$secondCousinTriple,
        alleleFreq=0.01), 0.002, tol=1e-3)
    expect_equal(RVsharing2(samplePedigrees$secondCousinPairWithLoop,
        alleleFreq=0.01), 0.058, tol=1e-3)
    expect_equal(RVsharing2(samplePedigrees$firstCousinInbreeding,
        alleleFreq=0.01), 0.122, tol=1e-3)
    expect_equal(RVsharing2(samplePedigrees$twoGenerationsInbreeding,
        alleleFreq=0.01), 0.029, tol=1e-3)
})

monteCarloComp <- function(ped, f)
{
    abs(RVsharing2(ped, alleleFreq=f) - 
        RVsharing2(ped, alleleFreq=f, nSimulations=1e4))
}

test_that('monte carlo sharing probabilities',
{
    data(samplePedigrees)
    tol <- 0.02
    set.seed(2)

    expect_equal(monteCarloComp(samplePedigrees$firstCousinPair),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPair),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$thirdCousinPair),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinTriple),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinTriple),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPairWithLoop),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinInbreeding),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$twoGenerationsInbreeding),
        0, tolerance=tol)

    expect_equal(monteCarloComp(samplePedigrees$firstCousinPair,
        0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPair,
        0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$thirdCousinPair,
        0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinTriple,
        0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinTriple,
        0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPairWithLoop,
        0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinInbreeding,
        0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$twoGenerationsInbreeding,
        0.01), 0, tolerance=tol)
})

