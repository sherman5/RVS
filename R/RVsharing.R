#' @include pedigree-methods.R
NULL

#' probability of sharing a rare variant among relatives
#' @export
#' @docType methods
#' @rdname RVsharing-methods
#'
#' @description computing probability that a rare variant is shared by a
#'  set of subjects in a pedigree using the gRain package
#' @details the function RVsharing computes the probability that all subjects
#'  identified as carriers of a rare variant in the vector carriers
#'  (or all final descendants in the pedigree if carriers == NULL) share that
#'  rare variant AND the final descendants not included in carriers do not
#'  carry it, given that the rare variant has been detected in any subject
#'  in the union of the carriers and the final descendants of the pedigree.
#'  A final descendant is defined as a subject without descendant in the
#'  pedigree, it it not necessarily in the youngest generation. If carriers
#'  enumerates a subset of pedigree members, the function will then compute
#'  the probability these carriers share the rare variant AND the final
#'  descendants not included in carriers do not carry it based on the above
#'  terms. To obtain the probability that a set of pedigree members carry a
#'  rare variant given it was seen in any of the set members (ignoring the
#'  carrier status of final descendants not in the set), the pedigree must be
#'  trimmed of the other final descendants before calling RVsharing.
#'
#' @param ped S3 pedigree object or a list of pedigree objects
#' @param carriers subjects in pedigree that have the variant, if
#'  ped is a list, then this will also be a list of vectors specifying
#'  the carriers in each pedigree
#' @param alleleFreq allele frequency among the founders
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @param nSim number of simulations used in monte carlo calculation
#' @param founderDist custom distribution among founders - only used
#'  when simulating probability with nSim
#' @param useAffected allows the user to condition on seeing the variant
#'  among the affected subjects instead of the final descendants
#' @param kinshipOrder order of the polynomial approximation to the distribtion
#'  of the number of distinct alleles in the founders (d in Bureau et al.).
#'  Must be <= 5
#' @param ... allows for arguments in the style of v1.7
#' @return sharing probability between all carriers in pedigree
#' @examples
#'  data("samplePedigrees")
#'  RVsharing(samplePedigrees$firstCousinPair)
#' @references Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E.,
#'  Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H.
#'  and Ruczinski, I. (2014) Inferring rare disease risk variants based on
#'  exact probabilities of sharing by multiple affected relatives.
#'  Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198.
setGeneric('RVsharing', function(ped, carriers=NULL, alleleFreq=NA,
kinshipCoeff=NA, nSim=NA, founderDist=NULL, useAffected=FALSE,
kinshipOrder=5, ...)
    {standardGeneric('RVsharing')})

#' @rdname RVsharing-methods
#' @aliases RVsharing
setMethod('RVsharing', signature(ped='pedigree'),
function(ped, carriers, alleleFreq, kinshipCoeff, nSim,
founderDist, useAffected, kinshipOrder, ...)
{
    # needed for backwards compatibility with v1.7
    ped <- oldArgs(ped, list(...)$data, list(...)$dad.id, list(...)$mom.id)

    # pre-processing step
    checkArgs(alleleFreq, kinshipCoeff, nSim, founderDist)
    if (!useAffected) ped$affected <- numeric(0)
    procPed <- processPedigree(ped, carriers)

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
        prob <- twoFounderSharingProb(procPed, kinshipCoeff, kinshipOrder)
    }
    else
    {
        prob <- oneFounderSharingProb(procPed)
    }
    
    # print and return result
    carrierText <- paste(procPed$origID[procPed$carriers], collapse=' ')
    affectedText <- paste(procPed$origID[procPed$affected], collapse=' ')
    message(paste('Probability subjects', carrierText, 'among',
        affectedText, 'share a rare variant:', signif(prob, 4)))
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
#'  if they are not
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
#'  the user called RVsharing using a function signature from v1.7, this
#'  will convert the arguments into a pedigree suitable for the
#'  signature in version > 2.0
#' @param ped a pedigree object
#' @param data numeric/character vector of subject ids
#' @param dad.id numeric/character vector of father ids, founders' parents
#'  should be NA or 0
#' @param mom.id numeric/character vector of mother ids, founders' parents
#'  should be NA or 0
#' @return if old arguments are provided, a pedigree object is returned,
#'  otherwise ped is returned
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
