context('multiple family calculations')

test_that('multipleFamilyPValue',
{
    # test specific result
    data(samplePedigrees)
    probs <- sapply(samplePedigrees, RVsharing)
    obs <- c(rep(TRUE, 2), rep(FALSE, 7))
    names(obs) <- names(probs)
    expect_equal(multipleFamilyPValue(probs, obs), 0.0135, tol=0.001)

    # test error handling
    expect_error(multipleFamilyPValue(unname(probs), obs))
    expect_error(multipleFamilyPValue(probs, unname(obs)))
    names(obs) <- LETTERS[1:length(probs)]
    expect_error(multipleFamilyPValue(probs, obs))
        
    # test agreement between each backend
    for (n in 1:100)
    {
        obs <- sample(c(TRUE, FALSE), length(probs), replace=TRUE)
        names(obs) <- names(probs)
        prob_R <- multipleFamilyPValue(probs, obs, backend='r')
        prob_CPP <- multipleFamilyPValue(probs, obs, backend='cpp')
        expect_equal(prob_R, prob_CPP)
    }
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

    # test consistency of backends across all variants
    result_R <- multipleVariantPValue(snpMat$genotypes, snpMat$fam, sharingProbs, backend='r')
    result_CPP <- multipleVariantPValue(snpMat$genotypes, snpMat$fam, sharingProbs, backend='cpp')
    expect_true(all(result_R$pvalues == result_CPP$pvalues))
    expect_true(all(result_R$potential_pvalues == result_CPP$potential_pvalues))
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

    # test specific result for both backends
    pval <- enrichmentPValue(snpMat$genotypes, snpMat$fam, sharingProbs, backend='r')
    expect_equal(pval, 0.124, tolerance=1e-3)
    pval <- enrichmentPValue(snpMat$genotypes, snpMat$fam, sharingProbs, backend='cpp')
    expect_equal(pval, 0.124, tolerance=1e-3)
})  

