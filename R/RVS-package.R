#' @keywords internal
#' @details
#' Rare Variant Sharing (RVS) implements tests of association and linkage between rare genetic variant genotypes and a dichotomous phenotype, e.g. a disease status, in family samples. The tests are based on probabilities of rare variant sharing by relatives under the null hypothesis of absence of linkage and association between the rare variants and the phenotype and apply to single variants or multiple variants in a region (e.g. gene-based test).
#' @importFrom stats runif weighted.mean
#' @importFrom utils combn data installed.packages
#' @import kinship2
#' @importClassesFrom snpStats SnpMatrix
#' @import gRain
#' @import methods
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
