context('Variable Elimination')

test_that('test creating individual pmf\'s',
{
    data(testPed)

    fam <- testPedList[[1]]
    fam$id <- 1:length(fam$id)

    margVars <- matrix(0, nrow=2,ncol=2)
    margVars[1,] <- c(1,9)
    pmf <- getIndividualPmf(3, c(4,5), margVars)

})
