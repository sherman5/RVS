#' @include pedigree-methods.R
NULL

#' probability of sharing a rare variant among relatives
#' @export
#' @docType methods
#' @rdname RVsharing-methods
#'
#' @description computing probability that a rare variant is shared by a
#' set of subjects in a pedigree using the gRain package
#' @details the function RVsharing computes the probability that all subjects
#' identified as carriers of a rare variant in the vector carriers
#' (or all final descendants in the pedigree if carriers == NULL) share that
#' rare variant AND the final descendants not included in carriers do not
#' carry it, given that the rare variant has been detected in any subject
#' in the union of the carriers and the final descendants of the pedigree,
#' when distinguishHomo = FALSE by default.
#' A final descendant is defined as a subject without descendant in the
#' pedigree, it it not necessarily in the youngest generation. If carriers
#' enumerates a subset of pedigree members, the function will then compute
#' the probability these carriers share the rare variant AND the final
#' descendants not included in carriers do not carry it based on the above
#' terms. To control the subset of pedigree members for which sharing of the
#' rare variant is evaluated, define the members of the subset as affected in
#' ped and set the option useAffected=TRUE. The only restriction is that all
#' carriers must be defined as affected. The probability of the (2 to the power
#' ncarriers) configurations of one or two copies of the variant allele among 
#' the carriers is instead computed when distinguishHomo = TRUE.
#' @param ped S3 pedigree object or a list of pedigree objects
#' @param carriers subjects in pedigree that have the variant, if
#' ped is a list, then this will also be a list of vectors specifying
#' the carriers in each pedigree
#' @param alleleFreq allele frequency among the founders
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @param nSim number of simulations used in monte carlo calculation
#' @param founderDist custom distribution among founders. Only used
#' when simulating probability with nSim
#' @param useAffected a logical value indicating whether to condition on seeing the variant
#' among the affected subjects in ped instead of the final descendants
#' @param kinshipOrder order of the polynomial approximation to the distribtion
#' of the number of distinct alleles in the founders (d in Bureau et al.).
#' Must be <= 5
#' @param splitPed a logical value indicating whether to split the pedigree in subpedigrees below each founder to enable computations in pedigrees too large to be stored in a single Bayesian network. Cannot be used in conjunction with alleleFreq or kinshipCoef. If useAffected=TRUE, this option produces correct probabilities only if all the affected subjects are non-founders of the pedigree. Erroneous probabilities may be returned if any founder is affected.
#' @param useFounderCouples a logical value indicating whether to exploit the interchangeability of the mother and father from founder couples to save computations. Warning! This works only when all founders have only one spouse. Set to FALSE if at least one founder has two or more spouses. Only used when splitPed = TRUE
#' @param distinguishHomo a logical value indicating whether to compute distinct probabilities for homozygous and heterozygous variant carrier status
#' @param ... allows for additional arguments
#' @return sharing probability between all carriers in pedigree, except when splitPed = TRUE, where a vector of sharing probabilities for all subsets of the carriers is returned, and when distinguishHomo = TRUE, where a vector of sharing probabilities of the (2 to the power ncarriers) configurations of one or two copies of the variant allele among the carriers is returned. In this case, the probabilities are named by binary strings, such that a "0" means 1 copy of the variant allele and a "1" means 2 copies of the variant allele.
#' @examples
#' data("samplePedigrees")
#' RVsharing(samplePedigrees$firstCousinPair)
#' @references Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E.,
#' Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H.
#' and Ruczinski, I. (2014) Inferring rare disease risk variants based on
#' exact probabilities of sharing by multiple affected relatives.
#' Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198.
#'
#' Sherman, T., Fu, J., Scharpf, R., Bureau, A., and Ruczinski, I. (2018)
#' Detection of rare disease variants in extended pedigrees using RVS.
#' Bioinformatics, 1-3, doi: 10.1093/bioinformatics/bty976
setGeneric('RVsharing', function(ped, carriers=NULL, alleleFreq=NA,
kinshipCoeff=NA, nSim=NA, founderDist=NULL, useAffected=FALSE,
kinshipOrder=5, splitPed=FALSE, useFounderCouples=TRUE, ncores=1, distinguishHomo=FALSE, ...)
    {standardGeneric('RVsharing')})

#' @rdname RVsharing-methods
#' @aliases RVsharing
setMethod('RVsharing', signature(ped='pedigree'),
function(ped, carriers, alleleFreq, kinshipCoeff, nSim,
founderDist, useAffected, kinshipOrder, splitPed, useFounderCouples, ncores, distinguishHomo, ...)
{
    # needed for backwards compatibility with v1.7
    ped <- oldArgs(ped, list(...)$data, list(...)$dad.id, list(...)$mom.id)

    # pre-processing step
    checkArgs(alleleFreq, kinshipCoeff, nSim, founderDist)
    if (!useAffected) ped$affected <- numeric(0)
    procPed <- processPedigree(ped, carriers)
    if (length(procPed$affected) < 2)
        stop('need at least 2 affected subjects')

    # calculate sharing prob with appropiate method
    if (!is.na(nSim))
    {
        prob <- monteCarloSharingProb(procPed=procPed, alleleFreq=alleleFreq,
            kinshipCoeff=kinshipCoeff, nSim=nSim, founderDist=founderDist,
            kinshipOrder=kinshipOrder)
    }
    else if (!is.na(alleleFreq))
    {
        prob <- exactSharingProb(procPed, alleleFreq)
    }        
    else if (!is.na(kinshipCoeff))
    {
        prob <- twoFounderSharingProb(procPed, kinshipCoeff, kinshipOrder, distinguishHomo)
        if (distinguishHomo) names(prob) = R.utils::intToBin(0:(2^length(procPed$carriers)-1))
    }
    else if (splitPed)
    {
        prob <- oneFounderSharingProbSplitting(procPed, useFounderCouples,ncores=ncores)
    }
    else
    {
        prob <- oneFounderSharingProb(procPed, distinguishHomo,ncores=ncores)
        if (distinguishHomo) names(prob) = R.utils::intToBin(0:(2^length(procPed$carriers)-1))
    }
    
    # print and return result
    carrierText <- paste(procPed$origID[procPed$carriers], collapse=' ')
    affectedText <- paste(procPed$origID[procPed$affected], collapse=' ')

    if (splitPed)
    {
        message('Probability every subset of subjects among ', affectedText,
            ' share a rare variant:')
        print(signif(prob, 4))
    }
    else
    {
        message(paste('Probability subjects', carrierText, 'among',
            affectedText, 'share a rare variant:', signif(prob[1], 4)))
    }
    return(prob)
})

#' @rdname RVsharing-methods
#' @aliases RVsharing
setMethod('RVsharing', signature(ped='list'),
function(ped, carriers, alleleFreq, kinshipCoeff, nSim,
founderDist, useAffected, kinshipOrder, ...)
{
    if (is.null(carriers)) carriers <- rep(NULL, length(ped))
    probs <- sapply(1:length(ped), function(i) RVsharing(ped[[i]],
        carriers[[i]], alleleFreq, kinshipCoeff, nSim, founderDist,
        kinshipOrder, ...))
    id <- as.character(sapply(ped, function(p) p$famid[1]))
    names(probs) <- id
    return(probs)
})

#' check arguments provided to RVsharing for validty
#' @keywords internal
#'
#' @description verifies that arguments are valid, throws an error
#' if they are not
#' @inheritParams RVsharing
#' @return throws error if arguments invalid
checkArgs <- function(alleleFreq, kinshipCoeff, nSim, founderDist)
{
    if (!is.na(alleleFreq) & !is.null(founderDist))
        warning('founderDist ignored since alleleFreq was provided')
    if (!is.na(kinshipCoeff) & !is.null(founderDist))
        warning('founderDist ignored since kinshipCoeff was provided')
    if (!is.na(alleleFreq) & !is.na(kinshipCoeff))
        stop('can\'t use both alleleFreq and kinshipCoeff')   
}

#' check for arguments in v1.7 format
#' @keywords internal
#'
#' @description check arguments provided in ... to see if
#' the user called RVsharing using a function signature from v1.7, this
#' will convert the arguments into a pedigree suitable for the
#' signature in version > 2.0
#' @param ped a pedigree object
#' @param data numeric/character vector of subject ids
#' @param dad.id numeric/character vector of father ids, founders' parents
#' should be NA or 0
#' @param mom.id numeric/character vector of mother ids, founders' parents
#' should be NA or 0
#' @return if old arguments are provided, a pedigree object is returned,
#' otherwise ped is returned
oldArgs <- function(ped, data, dad.id, mom.id)
{
    if (!is.null(data) & !is.null(dad.id) & !is.null(mom.id))
    {
        if (!is.element('kinship2', installed.packages()[,1]))
            stop('kinship2 package required when using dad.id/mom.id')
        data <- as.numeric(data)
        dad.id <- as.numeric(dad.id)
        dad.id[dad.id == NA] <- 0
        mom.id <- as.numeric(mom.id)
        mom.id[mom.id == NA] <- 0
        return(kinship2::pedigree(id=data, dadid=dad.id, momid=mom.id))
    }
    else if (is.null(data) & is.null(dad.id) & is.null(mom.id))
    {
        return(ped)
    }
    else
    {
        stop('need all 3 args: data, dad.id, mom.id')
    }
}
