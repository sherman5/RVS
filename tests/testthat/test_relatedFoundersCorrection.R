context('related founders correction')

test_that('inferNumAlleles',
{
    expect_equal(inferNumAlleles(0.01, 10), 19)
})

test_that('computePhiVec',
{
    expect_equal(computePhiVec(10), c(0.01085, 0.00541), tol=1e-3)
})

test_that('inferTheta',
{
    expect_equal(inferTheta(0.05, c(0.1, 0.1)), c(-1+sqrt(3), -1-sqrt(3)))
})

test_that('computePFU',
{
    expect_equal(computePFU(10, 0.01), 0.099, tol=0.01)
})

test_that('relatedFounderCorrection',
{
    expect_equal(relatedFoundersCorrection(10, 0.01), 0.893, tol=0.01)
})
