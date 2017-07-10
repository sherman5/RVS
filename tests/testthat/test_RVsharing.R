context('RVsharing top level')

test_that('args processed correctly',
{

})

test_that('old args processed correctly',
{


})

test_that('list of pedigrees',
{
    data(samplePedigrees)
    result <- RVsharing(samplePedigrees)
    print(result)

})

test_that('monte carlo close to exact',
{
    monteCarloComp <- function(ped, freq, kin)
    {
        abs(RVsharing(ped, alleleFreq=freq, kinshipCoeff=kin) - 
            RVsharing(ped, alleleFreq=freq, kinshipCoeff=kin,
            nSim=15000))
    }

    data(samplePedigrees)
    tol <- 0.01
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
        freq=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPair,
        freq=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$thirdCousinPair,
        freq=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinTriple,
        freq=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinTriple,
        freq=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPairWithLoop,
        freq=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinInbreeding,
        freq=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$twoGenerationsInbreeding,
        freq=0.01), 0, tolerance=tol)

    expect_equal(monteCarloComp(samplePedigrees$firstCousinPair,
        kin=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPair,
        kin=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$thirdCousinPair,
        kin=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinTriple,
        kin=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinTriple,
        kin=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPairWithLoop,
        kin=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinInbreeding,
        kin=0.01), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$twoGenerationsInbreeding,
        kin=0.01), 0, tolerance=tol)
})
