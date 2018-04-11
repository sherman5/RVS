context('multiple family calculations')

test_that('multiple family p-value',
{
    data(samplePedigrees)
    probs <- sapply(samplePedigrees, RVsharing)
    obs <- c(T,T,F,F,F,F,F,F,F)
    names(obs) <- names(probs)
    expect_equal(multipleFamilyPValue(probs, obs), 0.0135, tol=0.001)
})
