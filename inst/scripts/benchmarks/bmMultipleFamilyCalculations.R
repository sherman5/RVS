library(RVS)
library(snpStats)
library(hexbin)
data(for.exercise)
data(samplePedigrees)

timeToSeconds <- function(time)
{
    time$yday * 86400 + time$hour * 3600 + time$min * 60 + time$sec
}

benchmarkEnrichmentPValue <- function(snpMat, famInfo, sharingProbs, backend)
{
    start_time <- timeToSeconds(as.POSIXlt(Sys.time()))
    result <- RVS::enrichmentPValue(snpMat=snpMat,
                                    famInfo=famInfo,
                                    sharingProbs=sharingProbs,
                                    backend=backend)
    runtime <- timeToSeconds(as.POSIXlt(Sys.time())) - start_time
    message(paste0("Time (sec): ", round(runtime, 2)))
    print(result)
}

benchmarkFamilyPValue <- function(sharingProbs, backend='r')
{
    set.seed(123)
    start_time <- timeToSeconds(as.POSIXlt(Sys.time()))
    result <- c()
    for (n in 1:1000)
    {
        pattern <- sample(c(TRUE, FALSE), length(sharingProbs), replace=TRUE)
        names(pattern) <- names(sharingProbs)
        result <- c(result, RVS::multipleFamilyPValue(sharingProbs=sharingProbs,
                                            observedSharing=pattern,
                                            backend=backend))
    }
    runtime <- timeToSeconds(as.POSIXlt(Sys.time())) - start_time
    message(paste0("Time (sec): ", round(runtime, 2)))
    print(summary(result))
}

benchmarkVariantPValue <- function(snpMat, famInfo, sharingProbs, filter=NULL, alpha=0, backend='R')
{
    start_time <- timeToSeconds(as.POSIXlt(Sys.time()))
    result <- RVS::multipleVariantPValue(snpMat=snpMat,
                                    famInfo=famInfo,
                                    sharingProbs=sharingProbs,
                                    filter=filter,
                                    alpha=alpha,
                                    backend=backend)
    runtime <- timeToSeconds(as.POSIXlt(Sys.time())) - start_time
    message(paste0("Time (sec): ", round(runtime, 2)))
    print(summary(result))
    print(summary(result$pvalues))
    print(summary(result$potential_pvalues))
}

message("Creating Families")
A_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinPair)
B_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinPair)
C_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinTriple)
D_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinTriple)
fams <- c(A_fams, B_fams, C_fams, D_fams)
fam_names <- paste("fam", c('A', 'B', 'C', 'D'), sep="_")
famids <- apply(expand.grid(fam_names, 1:3), 1, function(row) paste0(row[1], row[2]))
for (i in 1:12)
{
    fams[[i]]$famid <- rep(famids[i], length(fams[[i]]$id))
}
sample <- snps.10[1:12,]
famInfo <- data.frame(pedigree=famids, affected=rep(2, 12))
sharingProbs <- suppressMessages(RVsharing(fams))

#message("Benchmarking multipleFamilyPValue, R Version")
#benchmarkFamilyPValue(sharingProbs, 'r')

#message("Benchmarking multipleFamilyPValue, C++ Version")
#benchmarkFamilyPValue(sharingProbs, 'cpp')

#message("Benchmarking enrichmentPValue, R Version")
#benchmarkEnrichmentPValue(sample, famInfo, sharingProbs, 'r')

#message("Benchmarking enrichmentPValue, C++ Version")
#benchmarkEnrichmentPValue(sample, famInfo, sharingProbs, 'cpp')

#message("Benchmarking multipleVariantPValue with default parameters, R Version")
#benchmarkVariantPValue(sample, famInfo, sharingProbs, backend='r')

message("Benchmarking multipleVariantPValue with default parameters, C++ Version")
benchmarkVariantPValue(sample, famInfo, sharingProbs, backend='cpp')

#message("Benchmarking multipleVariantPValue with bonferroni cutoff of alpha=0.05, R version")
#benchmarkVariantPValue(sample, famInfo, sharingProbs, "bonferroni", 0.05, 'r')

#message("Benchmarking multipleVariantPValue with bonferroni cutoff of alpha=0.05, C++ version")
#benchmarkVariantPValue(sample, famInfo, sharingProbs, "bonferroni", 0.05, 'cpp')

