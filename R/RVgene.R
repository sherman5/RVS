#' convert a list of SnpMatrices to a single matrix in a similiar format
#' as LINKAGE except with minor allele counts
#' @export
#'
#' @description creates a matrix in LINKAGE format using pedigree information
#' from a list of pedigree objects and genotype information from a list of
#' SnpMatrices
#' @param matList list of SnpMatrices
#' @param pedList list of pedigrees
#' @return matrix in LINKAGE format
#' @examples
#' data(samplePedigrees)
#' data(snpMat)
#' ped <- samplePedigrees$secondCousinTriple
#' ex.ped.mat <- SnpMatrixToCount(list(snpMat), list(ped))
SnpMatrixToCount <- function(matList, pedList)
{
    if (length(matList) != length(pedList))
        stop('number of pedigrees and SnpMatrices do not match')

    # convert to {NA,0,1,2}
    matList <- lapply(matList, function(mat) as(mat, 'numeric'))

    matList <- lapply(1:length(matList), function(i)
    {
        mat <- matrix(NA, nrow=length(pedList[[i]]$id), ncol=6)
        mat[,1] <- pedList[[i]]$famid
        mat[,2] <- pedList[[i]]$id
        mat[,6] <- pedList[[i]]$affected + 1

        ndx <- match(as.numeric(rownames(matList[[i]])), mat[,2])
        if (!all(!is.na(ndx)))
            stop('SnpMatrix has subjects not in pedigree')        
        mat <- mat[ndx,]
        return(cbind(mat, matList[[i]]))   
    })

    mat <- matList[[1]]
    if (length(matList) > 1)    
    {
        for (i in 2:length(matList))
        {
            mat <- merge(mat, matList[[i]], all=TRUE)
        }
    }
    return(mat)
}

#' extract carriers of minor allele
#' @keywords internal
#' 
#' @description finds the carriers of the minor allele at a specified site
#' @param ped pedigree coded in a ped file with either two alleles per
#' variant ("alleles"), or a count of one allele ("count")
#' @param site site where to record carriers
#' @param fam ID of the family for which to extract carriers
#' @param type representation of allele count
#' @param minor.allele id of minor allele
#' @return carriers in ped
extract_carriers = function(ped,site,fam,type="alleles",minor.allele=2)
{
    if (!(type %in% c("alleles","count"))) stop ("Invalid type ",type)

    if (type == "alleles")    
    {    
        if (missing(fam))
        {
            if (length(unique(ped[,1]))>1)
                stop(paste("More than one family in ped data and",
                    "no family specified."))
            genodat = ped[,5:6+2*site]
            if (any(genodat==0))
                stop ("Alleles can't be labelled 0.")
            subid = ped[,2]
        }    
        else
        {
            genodat = ped[ped[,1]==fam&ped[,6]==2,5:6+2*site]
            subid = ped[ped[,1]==fam&ped[,6]==2,2]        
        }

        carriers.bool = apply(genodat,1,function(vec)
            ifelse(any(is.na(vec)),FALSE,any(vec==minor.allele)) )
    }
    else # type == "count"
    {
        if (missing(fam))
        {
            if (length(unique(ped[,1]))>1)
                stop(paste("More than one family in ped data and",
                    "no family specified."))
            allelecount = ped[ped[,6]==2,6+site]
            subid = ped[,2]
        }
        else
        {
            allelecount = ped[ped[,1]==fam&ped[,6]==2,6+site]
            subid = ped[ped[,1]==fam&ped[,6]==2,2]        
        }
        carriers.bool = ifelse(is.na(allelecount),FALSE,allelecount>0)
    }
    # Return list of carriers
    as.character(subid[carriers.bool])
}

#' Probability of sharing of rare variants in a family sample within a gene
#' @export
#'
#' @description Computing probability of sharing of rare variants in
#' a family sample within a genomic region such as a gene.
#' @details The function extracts the carriers of the minor allele at each
#' entry in sites in each family where it is present in ped.mat (or in the
#' families specified in fams if that argument is specified). It then
#' computes exact rare variant sharing probabilities in each family for
#' each variant by calling \code{RVsharing}. If multiple rare variants are seen
#' in the same family, the smallest sharing probability among all rare
#' variants is retained. The joint rare variant sharing probability over
#' all families is obtained as the product of the family-specific
#' probabilities. The p-value of the test allowing for sharing by a subset
#' of affected subjects over the rare variants in the genomic region is
#' then computed as the sum of the probabilities of the possible
#' combinations of sharing patterns among all families with a probability
#' less than or equal to the observed joint probability and a total number
#' of carriers greater than or equal to the sum of the number of carriers
#' in all families, using the values in \code{pattern.prob.list}, \code{nequiv.list} and
#' \code{N.list}. The families where all affected subjects share a rare variant
#' are determined by verifying if the length of the carrier vector equals
#' the maximum value of \code{N.list} for that family. The p-value of the test
#' requiring sharing by all affected subjects is computed by calling
#' \code{multipleFamilyPValue}.
#'
#' @param data A list of \code{SnpMatrix} objects corresponding to each pedigree
#' object in ped.listfams, or a data.frame or matrix encoding
#' the pedigree information and genotype data in the standard LINKAGE ped
#' format or the PLINK raw format with additive component only 
#' (see PLINK web site [1]). From the pedigree information, only the family ID 
#' in the first column, the subject ID in the second column and the affection 
#' status in the sixth column are used
#' (columns 3 to 5 are ignored). Also, family members without genotype data do
#' not need to appear in this object. The genotype of each variant can be
#' coded in two ways, each corresponding to a different value of the type
#' option: a minor allele count on one column with missing values coded NA,
#' (type="count") or the identity of the two alleles on two consecutive columns,
#' with missing values coded 0 corresponding to the standard
#' LINKAGE ped format (type="alleles"). If you provide a \code{SnpMatrix} object
#' then the genotype should be coded as the minor allele count + 1, i.e. 01 is the
#' homozygous genotype for the common allele.
#' @param ped.listfams a list of \code{pedigree} objects, one object for each
#' pedigree for which genotype data are included in \code{data}.
#' @param sites a vector of the column indices of the variant sites to
#' test in \code{data}. If the argument fams is provided, the variant sites
#' are tested in each corresponding family in the fams vector (a variant
#' present in multiple families must then be repeated for every families
#' where it appears).
#' @param fams an optional character vector of the names of families
#' in \code{data} and \code{ped.listfams} carrying the variants listed in the 
#' corresponding position in \code{sites}. If missing, the names of the families 
#' carrying the minor allele at each position in \code{sites} are extracted from 
#' \code{data}
#' @param pattern.prob.list a list of precomputed rare variant sharing 
#' probabilities for all possible sharing patterns in the families in 
#' \code{data} and \code{ped.listfams} 
#' @param nequiv.list an optional vector of the number of configurations 
#' of rare variant sharing by the affected subjects corresponding to the 
#' same pattern and probability in \code{pattern.prob.list}. Default is a vector 
#' of 1s
#' @param N.list a vector of the number of affected subjects sharing a 
#' rare variant in the corresponding pattern in \code{pattern.prob.list}
#' @param type an optional character string taking value "alleles" or 
#' "count". Default is "alleles"
#' @param minor.allele.vec an optional vector of the minor alleles at each 
#' site in the \code{sites} vector. It is not needed if type="count". If it is 
#' missing and type="alleles", the minor allele is assumed to take the 
#' value 2
#' @param precomputed.prob an optional list of vectors precomputed rare 
#' variant sharing probabilities for families in \code{data} and \code{ped.listfams}. 
#' If the vectors are named, the names must be strings formed by the 
#' concatenation of the sorted carrier names separated by semi-columns. 
#' If the vectors are not 
#' named, the vectors must represent probabilities for all the possible 
#' values of \code{N.list} for the corresponding family (one probability per 
#' value of \code{N.list})
#' @param maxdim upper bound on the dimension of the array containing the 
#' joint distribution of the sharing patterns for all families in fams 
#' (to avoid running out of memory)
#' @param partial.sharing logical indicating whether the test allowing for sharing
#' by a subset of affected subjects should be performed. If FALSE, only 
#' the test requiring sharing
#' by all affected subjects is computed. Default is TRUE
#' @param ... other arguments to be passed to RVsharing
#' @return A list with items:
#' \code{p} P-value of the exact rare variant sharing test allowing for sharing
#' by a subset of affected subjects.
#' \code{pall} P-value of the exact rare variant sharing test requiring sharing
#' by all affected subjects.
#' \code{potentialp} Minimum achievable p-value if all affected subjects were
#' carriers of a rare variant.
#' @examples
#' data(samplePedigrees)
#' data(ex.ped.mat)
#' fam15157 <- samplePedigrees$secondCousinTriple
#' fam15157.pattern.prob = c(RVsharing(fam15157,carriers=c(15,16,17)),
#'     RVsharing(fam15157,carriers=c(15,16)),
#'     RVsharing(fam15157,carriers=c(15)))
#' fam15157.nequiv = c(1,3,3)
#' # check that distribution sums to 1
#' sum(fam15157.pattern.prob*fam15157.nequiv)
#' fam15157.N = 3:1
#' fam28003 <- samplePedigrees$firstAndSecondCousinsTriple
#' fam28003.pattern.prob = c(RVsharing(fam28003,carriers=c(36,104,110)),
#'     RVsharing(fam28003,carriers=c(36,104)),
#'     RVsharing(fam28003,carriers=c(104,110)),
#'     RVsharing(fam28003,carriers=c(36)),
#'     RVsharing(fam28003,carriers=c(104)))
#' fam28003.N = c(3,2,2,1,1)
#' fam28003.nequiv = c(1,2,1,1,2)
#' # check that distribution sums to 1
#' sum(fam28003.pattern.prob*fam28003.nequiv)
#' # Creating lists
#' ex.pattern.prob.list = list("15157"=fam15157.pattern.prob,"28003"=fam28003.pattern.prob)
#' ex.nequiv.list = list("15157"=fam15157.nequiv,"28003"=fam28003.nequiv)
#' ex.N.list = list("15157"=fam15157.N,"28003"=fam28003.N)
#' ex.ped.obj = list(fam15157,fam28003)
#' names(ex.ped.obj) = c("15157","28003")
#' sites = c(92,119)
#' minor.allele.vec=c(1,4)
#' RVgene(ex.ped.mat,ex.ped.obj,sites,
#'     pattern.prob.list=ex.pattern.prob.list,
#' nequiv.list=ex.nequiv.list,N.list=ex.N.list,
#'     minor.allele.vec=minor.allele.vec)
#' # calling with a SnpMatrix list
#' data(famVCF)
#' fam15157.snp = suppressWarnings(VariantAnnotation::genotypeToSnpMatrix(fam15157.vcf))
#' fam28003.snp = suppressWarnings(VariantAnnotation::genotypeToSnpMatrix(fam28003.vcf))
#' ex.SnpMatrix.list = list(fam15157=fam15157.snp$genotypes,fam28003=fam28003.snp$genotypes)
#' RVgene(ex.SnpMatrix.list,ex.ped.obj,sites,
#'     pattern.prob.list=ex.pattern.prob.list, nequiv.list=ex.nequiv.list,
#'     N.list=ex.N.list,minor.allele.vec=minor.allele.vec)
#' @references Bureau, A., Begum, F., Taub, M.A., Hetmanski, J., Parker, M.M.,
#' Albacha-Hejazi, H., Scott, A.F., et al. (2019) Inferring Disease Risk Genes
#' from Sequencing Data in Multiplex Pedigrees Through Sharing of Rare Variants.
#' Genet Epidemiol. 43(1):37-49. doi: 10.1002/gepi.22155.
RVgene <- function(data, ped.listfams, sites, fams, pattern.prob.list,
nequiv.list, N.list, type="alleles", minor.allele.vec,
precomputed.prob=list(0), maxdim = 1e9, partial.sharing=TRUE, ...)
{
    if (is(data, 'list'))
    {
        ped.mat <- SnpMatrixToCount(data, ped.listfams)
        type <- 'count'
    }
    else
    {
        ped.mat <- data
    }

    if (missing(nequiv.list))
    {
        nequiv.list = lapply(pattern.prob.list,function(vec) rep(1,length(vec)))
        names(nequiv.list) = names(pattern.prob.list)
    }
    
    if (type=="alleles")
    {    
        if (missing(minor.allele.vec)) minor.allele.vec = rep(2,length(sites))
        if (length(sites)!=length(minor.allele.vec))
            stop ("Lengths of sites and minor.allele.vec vectors differs.")
    }
    
    if (missing(fams))
    {
        fams.vec = sites.alongfams = NULL
        if (type=="alleles") 
        {
            minor.allele.alongfams = NULL
            for (i in 1:length(sites))
            {
                fams.site = unique(ped.mat[ped.mat[,6]==2 &
                    (ped.mat[,5+2*sites[i]]==minor.allele.vec[i] |
                    ped.mat[,6+2*sites[i]]==minor.allele.vec[i]),1])
                if (length(fams.site)==0) stop("No variant allele at site ",sites[i])
                if (is.factor(fams.site)) fams.site=as.character(fams.site)
                fams.vec = c(fams.vec,fams.site)
                sites.alongfams = c(sites.alongfams,
                    rep(sites[i],length(fams.site)))
                minor.allele.alongfams = c(minor.allele.alongfams,
                    rep(minor.allele.vec[i],length(fams.site)))
            }
        }
        else
        {
            for (i in 1:length(sites))
            {
                # Remove subjects with missing genotype
                ped.obs = ped.mat[!is.na(ped.mat[,6+sites[i]]),]
                fams.site = unique(ped.obs[ped.obs[,6]==2 &
                    ped.obs[,6+sites[i]]>0,1])
                if (length(fams.site)==0) stop("No variant allele at site ",sites[i])
                if (is.factor(fams.site)) fams.site=as.character(fams.site)
                fams.vec = c(fams.vec,fams.site)
                sites.alongfams = c(sites.alongfams,
                    rep(sites[i],length(fams.site)))            
            }
        }
    }
    else 
    {
        if (length(sites)!=length(fams))
            stop ("Lengths of fams and sites vectors differs.")
        fams.vec = fams
        sites.alongfams = sites
        if (type=="alleles") minor.allele.alongfams = minor.allele.vec
    }
            
    fams.vec = as.character(fams.vec)
    # Using the famid as name of the pedigree objects in ped.listfams
    fams.names = sapply(ped.listfams,function(fam)fam$famid[1])
    names(ped.listfams) = fams.names
    missing.fams = fams.vec[!(fams.vec%in%fams.names)] 
    if (length(missing.fams>0))
        stop ("Families ",missing.fams," not in ped.listfams.")
    missing.fams = fams.vec[!(fams.vec%in%names(pattern.prob.list))] 
    if (length(missing.fams>0))
        stop ("Families ",missing.fams," not in pattern.prob.list.")
    missing.fams = fams.vec[!(fams.vec%in%names(N.list))] 
    if (length(missing.fams>0))
        stop ("Families ",missing.fams," not in N.list.")
    
    famu = unique(fams.vec)
    famRVprob = famNcarriers = rep(NA,length(famu))
    names(famRVprob) = names(famNcarriers) = famu
    # Loop over the families
    for (f in 1:length(fams.vec))
    {
        # get carriers list
        if (type=="alleles")
        {
            carriers = extract_carriers(ped.mat,sites.alongfams[f],
                fams.vec[f], type="alleles",minor.allele.alongfams[f])
        }
        else
        {   
            carriers = extract_carriers(ped.mat,sites.alongfams[f],
                fams.vec[f],type=type)
        }
                
        # Computation of RV sharing probability
        if (length(carriers)>0) 
        {
            #cat (f,"\n")
            if (fams.vec[f] %in% names(precomputed.prob))
            {
                # If the precomputed probabilities for the current family
                # have no name, then assume the probabilities are listed
                # for each possible number of carriers in the family
                if (is.null(names(precomputed.prob[[
                fams.vec[f]]])))
                    tmp = precomputed.prob[[
                        fams.vec[f]]][length(carriers)]
                # Otherwise, the names are assumed to be carrier subsets 
                # separated by ; and the probability for the current carriers
                # is extracted
                else
                    tmp = precomputed.prob[[
                        fams.vec[f]]][paste(sort(carriers),collapse=";")]
            }
            else
                tmp = suppressMessages(RVsharing(ped.listfams[[fams.vec[f]]],
                    carriers=carriers,...))
            # If the RV has lower sharing probability, we keep it for this fam
            if (is.na(famRVprob[fams.vec[f]]) || (tmp > 0 & tmp < famRVprob[fams.vec[f]]))
            {
                famRVprob[fams.vec[f]] = tmp
                famNcarriers[fams.vec[f]] = length(carriers)
            }
        }
    }
    # Identify number of informative families
    fam.info = names(famRVprob)[!is.na(famRVprob)]
    nfam.info = length(fam.info)
    if (nfam.info>0) mdim = prod(sapply(N.list[fam.info],length))
    else mdim = 0

    if (partial.sharing)
    {
        if (mdim > maxdim) 
        {
            warning(paste("Number of possible combinations of sharing",
                "patterns is too high. Partial sharing test cannot be performed."))
            compute.p = FALSE       
        }
        else compute.p = TRUE
    }
    else  compute.p = FALSE
        
    # No informative family    
    if (nfam.info == 0) p = pall = potentialp = 1
    # One informative family
    else if (nfam.info == 1)
    {
        if (compute.p) p = sum((nequiv.list[[fam.info]]*pattern.prob.list[[fam.info]])
            [round(pattern.prob.list[[fam.info]],5) <= 
            round(famRVprob[fam.info],5) & N.list[[fam.info]] >= 
            famNcarriers[fam.info]])
        else p=NA
        potentialp = min(pattern.prob.list[[fam.info]])
        pall = ifelse(famNcarriers[fam.info]==max(N.list[[fam.info]]),
            potentialp,1)
    }
    else # > 1 informative family 
    {
    if (compute.p)
    {
        if (nfam.info == 2)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]])
            nequiv.array = outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]])
            N.array = outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+")
        }
        else if (nfam.info == 3)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]]),
                pattern.prob.list[[fam.info[3]]])
            nequiv.array = outer(outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]]),
                nequiv.list[[fam.info[3]]])
            N.array = outer(outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+"),
                N.list[[fam.info[3]]],"+")
        }
        else if (nfam.info == 4)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(outer(outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]]),
                pattern.prob.list[[fam.info[3]]]),
                pattern.prob.list[[fam.info[4]]])
            nequiv.array = outer(outer(outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]]),
                nequiv.list[[fam.info[3]]]),
                nequiv.list[[fam.info[4]]])    
            N.array = outer(outer(outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+"),
                N.list[[fam.info[3]]],"+"),
                N.list[[fam.info[4]]],"+")
        }
        else if (nfam.info == 5)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(outer(outer(outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]]),
                pattern.prob.list[[fam.info[3]]]),
                pattern.prob.list[[fam.info[4]]]),
                pattern.prob.list[[fam.info[5]]])
            nequiv.array = outer(outer(outer(outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]]),
                nequiv.list[[fam.info[3]]]),
                nequiv.list[[fam.info[4]]]),
                nequiv.list[[fam.info[5]]])    
            N.array = outer(outer(outer(outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+"),
                N.list[[fam.info[3]]],"+"),
                N.list[[fam.info[4]]],"+"),
                N.list[[fam.info[5]]],"+")
        } 
        else if (nfam.info == 6)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(outer(outer(outer(outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]]),
                pattern.prob.list[[fam.info[3]]]),
                pattern.prob.list[[fam.info[4]]]),
                pattern.prob.list[[fam.info[5]]]),
                pattern.prob.list[[fam.info[6]]])
            nequiv.array = outer(outer(outer(outer(outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]]),
                nequiv.list[[fam.info[3]]]),
                nequiv.list[[fam.info[4]]]),
                nequiv.list[[fam.info[5]]]),
                nequiv.list[[fam.info[6]]])    
            N.array = outer(outer(outer(outer(outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+"),
                N.list[[fam.info[3]]],"+"),
                N.list[[fam.info[4]]],"+"),
                N.list[[fam.info[5]]],"+"),
                N.list[[fam.info[6]]],"+")
        } 
        else if (nfam.info == 7)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(outer(outer(outer(outer(outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]]),
                pattern.prob.list[[fam.info[3]]]),
                pattern.prob.list[[fam.info[4]]]),
                pattern.prob.list[[fam.info[5]]]),
                pattern.prob.list[[fam.info[6]]]),
                pattern.prob.list[[fam.info[7]]])
            nequiv.array = outer(outer(outer(outer(outer(outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]]),
                nequiv.list[[fam.info[3]]]),
                nequiv.list[[fam.info[4]]]),
                nequiv.list[[fam.info[5]]]),
                nequiv.list[[fam.info[6]]]),
                nequiv.list[[fam.info[7]]])    
            N.array = outer(outer(outer(outer(outer(outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+"),
                N.list[[fam.info[3]]],"+"),
                N.list[[fam.info[4]]],"+"),
                N.list[[fam.info[5]]],"+"),
                N.list[[fam.info[6]]],"+"),
                N.list[[fam.info[7]]],"+")
        } 
        else if (nfam.info == 8)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(outer(outer(outer(outer(outer(outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]]),
                pattern.prob.list[[fam.info[3]]]),
                pattern.prob.list[[fam.info[4]]]),
                pattern.prob.list[[fam.info[5]]]),
                pattern.prob.list[[fam.info[6]]]),
                pattern.prob.list[[fam.info[7]]]),
                pattern.prob.list[[fam.info[8]]])
            nequiv.array = outer(outer(outer(outer(outer(outer(outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]]),
                nequiv.list[[fam.info[3]]]),
                nequiv.list[[fam.info[4]]]),
                nequiv.list[[fam.info[5]]]),
                nequiv.list[[fam.info[6]]]),
                nequiv.list[[fam.info[7]]]),
                nequiv.list[[fam.info[8]]])    
            N.array = outer(outer(outer(outer(outer(outer(outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+"),
                N.list[[fam.info[3]]],"+"),
                N.list[[fam.info[4]]],"+"),
                N.list[[fam.info[5]]],"+"),
                N.list[[fam.info[6]]],"+"),
                N.list[[fam.info[7]]],"+"),
                N.list[[fam.info[8]]],"+")
        } 
        else if (nfam.info == 9)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(outer(outer(outer(outer(outer(outer(
            outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]]),
                pattern.prob.list[[fam.info[3]]]),
                pattern.prob.list[[fam.info[4]]]),
                pattern.prob.list[[fam.info[5]]]),
                pattern.prob.list[[fam.info[6]]]),
                pattern.prob.list[[fam.info[7]]]),
                pattern.prob.list[[fam.info[8]]]),
                pattern.prob.list[[fam.info[9]]])
            nequiv.array = outer(outer(outer(outer(outer(outer(outer(outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]]),
                nequiv.list[[fam.info[3]]]),
                nequiv.list[[fam.info[4]]]),
                nequiv.list[[fam.info[5]]]),
                nequiv.list[[fam.info[6]]]),
                nequiv.list[[fam.info[7]]]),
                nequiv.list[[fam.info[8]]]),
                nequiv.list[[fam.info[9]]])     
            N.array = outer(outer(outer(outer(outer(outer(outer(outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+"),
                N.list[[fam.info[3]]],"+"),
                N.list[[fam.info[4]]],"+"),
                N.list[[fam.info[5]]],"+"),
                N.list[[fam.info[6]]],"+"),
                N.list[[fam.info[7]]],"+"),
                N.list[[fam.info[8]]],"+"),
                N.list[[fam.info[9]]],"+") 
        } 
        else if (nfam.info == 10)
        {
            # Creating matrices of joint probabilities, number of equivalent
            # patterns and number of carriers for the two informative families
            pattern.prob.array = outer(outer(outer(outer(outer(outer(outer(
            outer(outer(
                pattern.prob.list[[fam.info[1]]],
                pattern.prob.list[[fam.info[2]]]),
                pattern.prob.list[[fam.info[3]]]),
                pattern.prob.list[[fam.info[4]]]),
                pattern.prob.list[[fam.info[5]]]),
                pattern.prob.list[[fam.info[6]]]),
                pattern.prob.list[[fam.info[7]]]),
                pattern.prob.list[[fam.info[8]]]),
                pattern.prob.list[[fam.info[9]]]),
                pattern.prob.list[[fam.info[10]]])
            nequiv.array = outer(outer(outer(outer(outer(outer(outer(outer(
            outer(
                nequiv.list[[fam.info[1]]],
                nequiv.list[[fam.info[2]]]),
                nequiv.list[[fam.info[3]]]),
                nequiv.list[[fam.info[4]]]),
                nequiv.list[[fam.info[5]]]),
                nequiv.list[[fam.info[6]]]),
                nequiv.list[[fam.info[7]]]),
                nequiv.list[[fam.info[8]]]),
                nequiv.list[[fam.info[9]]]),
                nequiv.list[[fam.info[10]]])    
            N.array = outer(outer(outer(outer(outer(outer(outer(outer(outer(
                N.list[[fam.info[1]]],
                N.list[[fam.info[2]]],"+"),
                N.list[[fam.info[3]]],"+"),
                N.list[[fam.info[4]]],"+"),
                N.list[[fam.info[5]]],"+"),
                N.list[[fam.info[6]]],"+"),
                N.list[[fam.info[7]]],"+"),
                N.list[[fam.info[8]]],"+"),
                N.list[[fam.info[9]]],"+"),
                N.list[[fam.info[10]]],"+") 
        } 
        else 
        {
            warning ("More than 10 informative families.")
            compute.p = FALSE
        }
    }
    
        # Computing potential p-value
        potentialp = prod(sapply(pattern.prob.list[fam.info],function (vec) min(vec[vec>0])))
        # Computing p-value
        pobs =  round(prod(famRVprob[fam.info]),5)
        if (compute.p)
            p = sum((nequiv.array*pattern.prob.array)[round(pattern.prob.array,
                5)<=pobs & N.array>=sum(famNcarriers[fam.info])])
        else p = NA
        maxN = sapply(N.list[fam.info],max)
        not = fam.info[famNcarriers[fam.info]<maxN]
        if (length(not)>0)
        {
            if (2^nfam.info <= maxdim)
            {
                pshare = list(ped.tocompute.vec=fam.info,pshare=
                    sapply(pattern.prob.list[fam.info],function (vec) min(vec[vec>0])))
                pall = suppressWarnings(get.psubset(fam.info,not,pshare))
            }
            else pall = NA
        }
        else pall = potentialp
    }
    list(p=p,pall=pall,potentialp=potentialp)    
}
