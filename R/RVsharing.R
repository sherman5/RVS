#' @include pedigree-methods.R
NULL

#' compute rare variant sharing probabilities
#' @export
#' @docType methods
#' @rdname RVsharing-methods
#'
#' @description calculates rare variant sharing probability between
#'  carriers in a pedigree given that the variant in seen in one of the
#   affected subjects
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
#'  Important note: the affected element of the pedigree object is ignored by
#'  RVsharing.
#'
#' @param ped S3 pedigree object or a list of pedigree objects
#' @param carriers #TODO
#' @param alleleFreq allele frequency among the founders
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @param nSim number of simulations used in monte carlo calculation
#' @param founderDist #TODO
#' @param ... #TODO
#' @return sharing probability between all carriers in pedigree
#' @examples
#'  data("samplePedigrees")
#'  RVsharing(samplePedigrees[[1]])
setGeneric('RVsharing', function(ped, carriers, alleleFreq, kinshipCoeff,
nSim, founderDist, ...) {standardGeneric('RVsharing')})

#' @rdname RVsharing-methods
#' @aliases #TODO
setMethod('RVsharing', signature(ped='pedigree'),
function(ped, carriers, alleleFreq, kinshipCoeff, nSim, founderDist, ...)
{
    # needed for backwards compatibility with v1.7
    ped <- oldArgs(ped, list(...)$data, list(...)$dad.id, list(...)$mom.id)

    # pre-processing step
    checkArgs(alleleFreq, kinshipCoeff, nSim, founderDist)
    procPed <- processPedigree(ped, carriers)

    # calculate sharing prob with appropiate method
    if (!missing(nSim))
    {
        return(monteCarloSharingProb(procPed=procPed, alleleFreq=alleleFreq,
            kinshipCoeff=kinshipCoeff, nSim=nSim, founderDist=founderDist))
    }
    else if (!missing(alleleFreq))
    {
        return(exactSharingProb(procPed, alleleFreq))
    }        
    else if (!missing(kinshipCoeff))
    {
        return(twoFounderSharingProb(procPed, kinshipCoeff))
    }
    else
    {
        return(oneFounderSharingProb(procPed))
    }
})

#' @rdname RVsharing-methods
#' @aliases #TODO
setMethod('RVsharing', signature(ped='list'),
function(ped, carriers, alleleFreq, kinshipCoeff, nSim, founderDist, ...)
{
    sapply(1:length(ped), function(i)
    {
        prob <- RVsharing(ped[[i]], carriers[[i]], alleleFreq,
            kinshipCoeff, nSim, founderDist, ...)
        if (!is.null(ped[[i]]@famid)) id <- as.character(ped[[i]]@famid)
        else id <- 'no_id'
        return(c(id, prob))
    })
})

#' check arguments provided to RVsharing for validty
#' @keywords internal
#'
#' @description verifies that arguments are valid, throws an error
#'  if they are not
#' @param alleleFreq see RVsharing
#' @param kinshipCoeff see RVsharing
#' @param nSim see RVsharing
#' @param founderDist see RVsharing
#' @param ... see RVsharing
#' @return throws error if arguments invalid
checkArgs <- function(alleleFreq, kinshipCoeff, nSim, founderDist)
{
    if (!missing(alleleFreq) & !missing(founderDist))
        warning('founderDist ignored since alleleFreq was provided')
    if (!missing(alleleFreq) & !missing(kinshipCoeff))
        stop('can\'t use both alleleFreq and kinshipCoeff')   
}

#' check for arguments in v1.7 format
#' @keywords internal
#'
#' @description check arguments provided in ... to see if
#'  the user called RVsharing using a function signature from v1.7, this
#'  will convert the arguments into a pedigree suitable for the
#'  signature in version > 2.0
#' @param ped pedigree object
#' @param data old argument
#' @param dad.id old argument
#' @param mom.id old argument
#' @return if old arguments are provided, a pedigree object is returned,
#'  otherwise ped is returned
oldArgs <- function(ped, data, dad.id, mom.id)
{
    if (!is.null(data) & !is.null(dad.id) & !is.null(mom.id))
    {
        if (!is.element('kinship2', installed.packages()[,1]))
            stop('kinship2 package required when using dad.id/mom.id')
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
