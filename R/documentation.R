#' RVS
#'
#' @docType package
#' @name RVS
#' @description Rare Variant Sharing (RVS) implements tests of association and linkage between rare genetic variant genotypes and a dichotomous phenotype, e.g. a disease status, in family samples. The tests are based on probabilities of rare variant sharing by relatives under the null hypothesis of absence of linkage and association between the rare variants and the phenotype and apply to single variants or multiple variants in a region (e.g. gene-based test).
#' @importFrom stats runif weighted.mean
#' @importFrom utils combn data installed.packages
#' @import kinship2
#' @importClassesFrom snpStats SnpMatrix
#' @import gRain
#' @import methods
NULL

#' list of 8 sample pedigree objects
#' 
#' @docType data
#' @name samplePedigrees
#' @usage samplePedigrees
NULL

#' matrix of pedigree information and genotype data from famVCF stored in the LINKAGE format
#' 
#' @docType data
#' @name ex.ped.mat
#' @usage ex.ped.mat
NULL

#' SnpMatrix with genotype information from famVCF for fam15157
#' 
#' @docType data
#' @name snpMat
#' @usage snpMat
NULL

#' VCF objects containing genotype data for two families: fam15157 and fam28003 (corresponding to the secondCousinTriple and firstAndSecondCousinsTriple families in samplePedigrees)
#' 
#' @docType data
#' @name fam28003.vcf
#' @usage fam28003.vcf
NULL

#' VCF objects containing genotype data for two families: fam15157 and fam28003 (corresponding to the secondCousinTriple and firstAndSecondCousinsTriple families in samplePedigrees)
#' 
#' @docType data
#' @name fam15157.vcf
#' @usage fam15157.vcf
NULL
