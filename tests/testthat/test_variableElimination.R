context('Variable Elimination')

    data(testPed)

    fam <- testPedList[[1]]
    fam$id <- 1:length(fam$id)


test_that('test creating individual pmf\'s',
{
    margVars <- matrix(0, nrow=2,ncol=1)
    margVars[1,] <- c(4)

    pmf <- getIndividualPmf(1, c(2,3), margVars)
    print(pmf)
    expect_equal(pmf[['prob']][1,1,1], 1)
    expect_equal(pmf[['prob']][1,3,3], 0)
    expect_equal(pmf[['prob']][1,2,2], 0.25)
    expect_equal(pmf[['prob']][2,2,2], 0.5)
    expect_equal(pmf[['prob']][3,2,2], 0.25)
})

test_that('test creating individual pmf\'s - plugging in marginal variables',
{

})
