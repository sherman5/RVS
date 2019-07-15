context('multiple family calculations')

test_that('multiple family p-value',
{
    data(samplePedigrees)
    probs <- sapply(samplePedigrees, RVsharing)
    obs <- c(rep(TRUE, 2), rep(FALSE, 7))
    names(obs) <- names(probs)
    expect_equal(multipleFamilyPValue(probs, obs), 0.0135, tol=0.001)
})


test_that('multiple variant p-value',
{
    pedfile <- system.file("extdata/sample.ped.gz", package="RVS")
    sample <- snpStats::read.pedfile(pedfile, snps=paste('variant', LETTERS[1:20], sep='_'))
    A_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinPair)
    B_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinPair)
    C_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinTriple)
    D_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinTriple)
    fams <- c(A_fams, B_fams, C_fams, D_fams)
    famids <- unique(sample$fam$pedigree)
    for (i in 1:12)
    {
        fams[[i]]$famid <- rep(famids[i], length(fams[[i]]$id))
    }
    sharingProbs <- suppressMessages(RVsharing(fams))
    result <- multipleVariantPValue(sample$genotypes, sample$fam, sharingProbs)
    expect_equal(length(result$pvalues), 20)
})

