library(RVS)

timeToSeconds <- function(time)
{
    time$yday * 86400 + time$hour * 3600 + time$min * 60 + time$sec
}

benchmarkEnrichmentPValue <- function(snpMat, famInfo, sharingProbs, backend)
{
    start_time <- timeToSeconds(as.POSIXlt(Sys.time()))
    result <- enrichmentPValue(snpMat=snpMat,
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
        result <- c(result, multipleFamilyPValue(sharingProbs=sharingProbs,
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
    result <- multipleVariantPValue(snpMat=snpMat,
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

# create family data and calculate sharing probs
data(samplePedigrees)
A_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinPair)
B_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinPair)
C_fams <- lapply(1:3, function(i) samplePedigrees$firstCousinTriple)
D_fams <- lapply(1:3, function(i) samplePedigrees$secondCousinTriple)
fams <- c(A_fams, B_fams, C_fams, D_fams)
fam_names <- paste("fam", c('A', 'B', 'C', 'D'), sep="_")
famids <- apply(expand.grid(fam_names, 1:3), 1, function(row) paste0(row[1], row[2]))
famInfo <- data.frame()
for (i in 1:12)
{
    fams[[i]]$famid <- rep(famids[i], length(fams[[i]]$id))
    famInfo <- rbind(famInfo, data.frame(
        pedigree=fams[[i]]$famid,
        affected=fams[[i]]$affected
    ))
}
famInfo <- famInfo[famInfo$affected==1,]
sharingProbs <- suppressMessages(RVsharing(fams))

library(snpStats)
data(for.exercise)
#snpMat <- snps.10[1:nrow(famInfo),4]
snpMat <- snps.10[1:nrow(famInfo),1:1000]

#pedfile <- system.file("extdata/sample.ped.gz", package="RVS")
#sample <- snpStats::read.pedfile(pedfile, snps=paste('variant', LETTERS[1:20], sep='_'))
#snpMat <- sample$genotypes[1:nrow(famInfo),]

#message("Benchmarking multipleFamilyPValue, R Version")
#benchmarkFamilyPValue(sharingProbs, 'r')

#message("Benchmarking multipleFamilyPValue, C++ Version")
#benchmarkFamilyPValue(sharingProbs, 'cpp')

#message("Benchmarking enrichmentPValue, R Version")
#benchmarkEnrichmentPValue(sample, famInfo, sharingProbs, 'r')

#message("Benchmarking enrichmentPValue, C++ Version")
#benchmarkEnrichmentPValue(sample, famInfo, sharingProbs, 'cpp')

message("Benchmarking multipleVariantPValue with default parameters, R Version")
benchmarkVariantPValue(snpMat, famInfo, sharingProbs, backend='r')

message("Benchmarking multipleVariantPValue with default parameters, C++ Version")
benchmarkVariantPValue(snpMat, famInfo, sharingProbs, backend='cpp')

result_R <- multipleVariantPValue(snpMat=snpMat, famInfo=famInfo, sharingProbs=sharingProbs,
                                  filter=NULL, alpha=0, backend='r')

result_CPP <- multipleVariantPValue(snpMat=snpMat, famInfo=famInfo, sharingProbs=sharingProbs,
                                    filter=NULL, alpha=0, backend='cpp')

#print(result_R)
#print(result_CPP)
print(all.equal(unname(result_R$pvalues), unname(result_CPP$pvalues)))

#message("Benchmarking multipleVariantPValue with bonferroni cutoff of alpha=0.05, R version")
#benchmarkVariantPValue(sample, famInfo, sharingProbs, "bonferroni", 0.05, 'r')

#message("Benchmarking multipleVariantPValue with bonferroni cutoff of alpha=0.05, C++ version")
#benchmarkVariantPValue(snpMat, famInfo, sharingProbs, "bonferroni", 0.05, 'cpp')

