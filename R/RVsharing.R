#' Compute Rare Variant Sharing Probabilities
#' @export
#'
#' @description Calculates rare variant sharing probability between
#'  carriers in a pedigree given that the variant in seen in one of the
#   affected subjects
#' @details
#'
#' @param ped S3 Pedigree object
#' @param alleleFreq allele frequency among the founders
#' @param kinshipCoeff mean kinship coefficient among the founders
#' @param nSim number of simulations used in monte carlo calculation
#' @return sharing probability between all carriers in pedigree
#' @examples
#'  data("samplePedigrees")
#'  RVsharing(samplePedigrees[[1]])
setGeneric('RVsharing', function(ped, carriers, alleleFreq, kinshipCoeff,
nSim, founderDist, ...) {standardGeneric('RVsharing')})

setMethod('RVsharing', signature(),
function(ped, carriers, alleleFreq, kinshipCoeff, nSim, founderDist, ...)
{
    # needed for backwards compatibility with v1.7
    ped <- oldArgs(ped, list(...)$data, list(...)$dad.id, list(...)$mom.id)

    # pre-processing step
    checkArgs(alleleFreq, kinshipCoeff, nSim, founderDist, ...)
    data(mendelProbTable)
    procPed <- processPedigree(ped, carriers)

    # calculate sharing prob with appropiate method
    if (!missing(nSim))
    {
        return(monteCarloSharingProb(procPed, alleleFreq, kinshipCoeff,
            nSim, founderDist))
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

setMethod('RVsharing', signature(ped='list'),
function(ped, carriers, alleleFreq, kinshipCoeff, nSim, founderDist, ...)
{
    sapply(1:length(ped), function(i) RVsharing(ped[[i]], carriers[[i]],    
        alleleFreq, kinshipCoeff, nSim, founderDist, ...))
})

checkArgs <- function(alleleFreq, kinshipCoeff, nSim, founderDist, ...)
{
    if (!missing(alleleFreq) & !missing(founderDist))
        warning('founderDist ignored since alleleFreq was provided')
    if (!missing(alleleFreq) & !missing(kinshipCoeff))
        stop('can\'t use both alleleFreq and kinshipCoeff')   
}

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
