context('multiple family calculations')

test_that('multiple family p-value',
{
    data(samplePedigrees)
    probs <- sapply(samplePedigrees, RVsharing)
    obs <- c(rep(TRUE, 2), rep(FALSE, 7))
    names(obs) <- names(probs)
    expect_equal(multipleFamilyPValue(probs, obs), 0.0135, tol=0.001)
})
