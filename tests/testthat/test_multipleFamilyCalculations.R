context('multiple family calculations')

test_that('multiple family p-value',
{
    data(samplePedigrees)
    probs <- sapply(samplePedigrees, RVsharing)
    expect_equal(multipleFamilyPValue(probs, c(T,T,F,F,F,F,F,F,F)), 0.00551,
        tol=0.001)
})
