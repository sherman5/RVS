#' RVS
#'
#' @docType package
#' @name RVS
#' @importFrom stats runif weighted.mean
#' @importFrom utils combn data installed.packages
#' @import kinship2
#' @importClassesFrom snpStats SnpMatrix
#' @import gRain
#' @import methods
#' @useDynLib RVS
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
