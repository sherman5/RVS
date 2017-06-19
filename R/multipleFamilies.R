#' Probability of sharing of rare variants in a subset of families
#' @export
#'
#' @description Computing probability of sharing of rare variants in
#'  a subset of families where rare variants are seen based on precomputed
#'  family-specific rare variant sharing probabilities.
#' @details All the subsets of families of size equal or inferior to the
#'  length of not are created, and the joint probability of each such
#'  subset not sharing a rare variant and the remaining families sharing
#'  a rare variant is obtained as the product of the family-specific rare
#'  variant sharing probabilities or its complement. The function then sums
#'  the pattern probabilities inferior or equal to the probability
#'  of the observed pattern of the not families not sharing a rare variant
#'  and the remaining families sharing a rare variant.
#' @param shared vector of sharing probabilities for families where
#'  all affected subjects share the variant
#' @param notShared vector of sharing probabilities for families where
#'  not all affected subjects share the variant
#' @return p-value
multipleFamilyPValue <- function(shared, notShared)
{
    # check: "not" contains at least one family 
    if (length(notShared)==0)
        stop("Vector 'notShared' of families not sharing the RV is empty.")

    # If all families share the variant, then return 1
#    if (length(not)==length(vec)) return (1)

    vec <- 1:(length(notShared)+length(shared))
    not <- 1:length(notShared)
    pshare.data <- data.frame(pshare=c(notShared,shared),
        ped.tocompute.vec=vec)
    
    p.vec = pshare.data$pshare[pshare.data$ped.tocompute.vec%in%vec]
    names(p.vec) = pshare.data$ped.tocompute.vec[pshare.data$ped.tocompute.vec%in%vec]
    nf = length(p.vec)
    nnot = sum(names(p.vec)%in%not)

    # Probability of observed data
    p.obs = prod(p.vec[!(names(p.vec)%in%not)],1-p.vec[as.character(not)])

    # Tail probability includes case where all families share the variant
    p = prod(p.vec)
    for (h in 1:nnot)
    {
        comb.mat = combn(nf,h)
        #	print(comb.mat)
        # plus the cases where the probability with two families not sharing the variant is less extreme than the observed
        # Compute probability for all pairs of families not sharing
        for (i in 1:ncol(comb.mat))
        {
            ptmp = prod(p.vec[-comb.mat[,i]],1-p.vec[comb.mat[,i]])
            #    print(ptmp)
            if (ptmp <= p.obs) p = p + ptmp
        }
    }
    return(p)
}
