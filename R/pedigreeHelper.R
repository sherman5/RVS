#' \code{processPedigree} pre-process a pedigree for calculations
#' 
#' @param ped pedigree objevt (S3)
#' @return list containing relevant pedigree info
#' @keywords internal
processPedigree <- function(ped, carriers)
{
    carriers <- which(ped$id %in% carriers)
    ped$id <- 1:length(ped$id)
    parents <- sapply(ped$id, function(i) c(ped$findex[i], ped$mindex[i]))
    founders <- which(parents[1,] == 0)
    affected <- which(ped$affected == 1)

    if (sum(affected %in% founders) > 0)
        stop('some founders are affected')
    if (sum(ped$affected==1) < 2)
        stop('need at least 2 affected subjects')

    return(list('ped'=ped, 'parents'=parents, 'founders'=founders,
        'affected'=affected, 'size'=length(ped$id),
        'carriers'=carriers, 'id'=ped$id))
}

