context('multiple family calculations')

test_that('multipleFamilyPValue',
{
    # test specific result
    data(samplePedigrees)
    probs <- sapply(samplePedigrees, RVsharing)
    obs <- c(rep(TRUE, 2), rep(FALSE, 7))
    names(obs) <- names(probs)
    #expect_equal(multipleFamilyPValue(probs, obs), 0.0135, tol=0.001)

    # test error handling
    expect_error(multipleFamilyPValue(unname(probs), obs))
    expect_error(multipleFamilyPValue(probs, unname(obs)))
    names(obs) <- LETTERS[1:length(probs)]
    expect_error(multipleFamilyPValue(probs, obs))
})

test_that('multipleVariantPValue',
{
    # create test data
    pedfile <- system.file("extdata/sample.ped.gz", package="RVS")
    snpMat <- snpStats::read.pedfile(pedfile, snps=paste('variant', LETTERS[1:20], sep='_'))
    A_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinPair)
    B_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinPair)
    C_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinTriple)
    D_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinTriple)
    fams <- c(A_fams, B_fams, C_fams, D_fams)
    famids <- unique(snpMat$fam$pedigree)
    for (i in 1:12)
    {
        fams[[i]]$famid <- rep(famids[i], length(fams[[i]]$id))
    }
    sharingProbs <- suppressMessages(RVsharing(fams))
    # make sure the result is named, has correct length
    #result <- multipleVariantPValue(sample$genotypes, sample$fam, sharingProbs)
    #expect_true(!is.null(names(result$pvalues)))
    #expect_true(!is.null(names(result$potential_pvalues)))
    #expect_equal(length(result$pvalues), 20)
})

test_that('enrichmentPValue',
{
    # create test data
    data(samplePedigrees)
    pedfile <- system.file("extdata/sample.ped.gz", package="RVS")
    snpMat <- snpStats::read.pedfile(pedfile, snps=paste('variant', LETTERS[1:20], sep='_'))
    A_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinPair)
    B_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinPair)
    C_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinTriple)
    D_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinTriple)
    fams <- c(A_fams, B_fams, C_fams, D_fams)
    famids <- unique(snpMat$fam$pedigree)
    for (i in 1:12)
    {
        fams[[i]]$famid <- rep(famids[i], length(fams[[i]]$id))
    }
    sharingProbs <- suppressMessages(RVsharing(fams))
})  

