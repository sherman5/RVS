context('RVsharing top level')

test_that('args processed correctly',
{
    data(samplePedigrees)
    ped <- samplePedigrees[[1]]
    f <- function(n) rep(1,n)

    expect_error(RVsharing(ped, alleleFreq=0.1, kinshipCoeff=0.1))
    expect_warning(RVsharing(ped, alleleFreq=0.1, founderDist=f))
    expect_warning(RVsharing(ped, kinshipCoeff=0.1, founderDist=f))
})

test_that('list of pedigrees',
{
    data(samplePedigrees)

    probs <- sapply(samplePedigrees, RVsharing)
    ids <- as.character(sapply(samplePedigrees, function(p) p$famid[1]))
    result <- RVsharing(samplePedigrees)

    expect_equal(names(result), ids)
    expect_equal(unname(result), unname(probs))
})

test_that('monte carlo close to exact',
{
    set.seed(0)

    monteCarloComp <- function(ped, freq=NA, kin=NA)
    {
        abs(RVsharing(ped, alleleFreq=freq, kinshipCoeff=kin) - 
            RVsharing(ped, alleleFreq=freq, kinshipCoeff=kin, nSim=2e4))
    }

    data(samplePedigrees)
    tol <- 0.02
    MAF <- 0.01

    expect_equal(monteCarloComp(samplePedigrees$firstCousinPair),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPair),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPairWithLoop),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinInbreeding),
        0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$twoGenerationsInbreeding),
        0, tolerance=tol)

    expect_equal(monteCarloComp(samplePedigrees$firstCousinPair,
        freq=MAF), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinTriple,
        freq=MAF), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$secondCousinPairWithLoop,
        freq=MAF), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$firstCousinInbreeding,
        freq=MAF), 0, tolerance=tol)
    expect_equal(monteCarloComp(samplePedigrees$twoGenerationsInbreeding,
        freq=MAF), 0, tolerance=tol)
})
