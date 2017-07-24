#' probability of sharing of rare variants among all affected relatives in
#' a family sample within a gene
#' @export
#'
#' @description Computing probability of sharing of rare variants among all
#'  affected relatives in a family sample within a genomic region such as
#'  a gene.
#' @details The function extracts the carriers of the minor allele at each 
#'  entry in sites in each family where it is present in ped.mat (or in the 
#'  families specified in fams if that argument is specified). It then 
#'  computes exact rare variant sharing probabilities in each family for 
#'  each variant by calling RVsharing. If multiple rare variants are seen 
#'  in the same family, the smallest sharing probability among all rare 
#'  variants is retained. The families where all affected subjects share a 
#'  rare variant are determined by verifying if the length of the carrier 
#'  vector equals the number of affected subjects that family. The joint 
#'  rare variant sharing probability over all families is obtained as the 
#'  product of the family-specific probabilities. The p-value of the test 
#'  requiring sharing by all affected subjects is computed by calling 
#'  multipleFamilyPValue
#'
#' @param data A list of SnpMatrix objects corresponding to each pedigree
#'  object in ped.listfams, alternatively a data.frame or matrix encoding
#'  the pedigree information and genotype data in the standard LINKAGE ped
#'  format (see PLINK web site [1]). In fact, only the family ID in the first
#'  column, the subject ID in the second column, the affection status in the
#'  sixth column and the genotype data starting in the seventh column are used
#'  (columns 3 to 5 are ignored). Also, family members without genotype data do
#'  not need to appear in this matrix. The genotype of each variant can be
#'  coded in two ways, each corresponding to a different value of the type
#'  option: a minor allele count on one column, as returned for example by the
#'  genotypeToSnpMatrix function, with missing values coded NA
#'  (type="count") or the identity of the two alleles on two consecutive
#'  columns, with missing values coded 0 (type="alleles")
#' @param ped.listfams a list of pedigree objects, one object for each 
#'  pedigree in ped.mat
#' @param sites a vector of the column indices of the variant sites to test 
#'  in ped.mat. If the argument fams is provided, the variant sites are 
#'  tested in each corresponding family in the fams vector (a variant 
#'  present in multiple families must then be repeated for every families 
#'  where it appears)
#' @param fams an optional character vector of the names of families in 
#'  ped.mat and ped.listfams carrying the variants listed in the 
#'  corresponding position in sites. If missing, the names of the families 
#'  carrying the minor allele at each position in sites are extracted from 
#'  ped.mat.
#' @param pshare.vec a vector of the probabilities that all affected 
#'  relatives share a rare variant, for every family. This vector must be 
#'  named with the family names
#' @param type an optional character string taking value "alleles" or 
#'  "count". Default is "alleles"
#' @param minor.allele.vec an optional vector of the minor alleles at each 
#'  site in the sites vector. It is not needed if type="count". If it is 
#'  missing and type="alleles", the minor allele is assumed to take the 
#'  value 2
#' @param precomputed.prob an optional list of vectors precomputed rare 
#'  variant sharing probabilities for families in ped.mat and ped.listfams.
#'  The vectors represent probabilities for all the possible values of 
#'  N.list for the corresponding family (one probability per value of 
#'  N.list)
#' @param maxdim upper bound on the dimension of the array containing the 
#'  joint distribution of the sharing patterns for all families in fams 
#'  (to avoid running out of memory)
#' @return A list with items:
#'  pall P-value of the exact rare variant sharing test requiring sharing 
#'  by all affected subjects.
#'  potentialp Minimum achievable p-value if all affected subjects were 
#'  carriers of a rare variant.
#' @references http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped	
#' @references Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E.,
#'  Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H.
#'  and Ruczinski, I. (2014) Inferring rare disease risk variants based on
#'  exact probabilities of sharing by multiple affected relatives.
#'  Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198.
RVgene_allshare <- function(data, ped.listfams, sites, fams, pshare.vec,
type="alleles", minor.allele.vec, precomputed.prob=list(0), maxdim = 1e9)
{
    if (class(data) == 'list')
    {
        ped.mat <- SnpMatrixToCount(data, ped.listfams)
        type <- 'count'
    }
    else
    {
        ped.mat <- data
    }

    if (type=="alleles")
    {    
        if (missing(minor.allele.vec)) minor.allele.vec = rep(2,length(sites))    
        if (length(sites)!=length(minor.allele.vec)) stop ("Lengths of sites and minor.allele.vec vectors differs.")
    }
    
    if (missing(fams))
    {
        fams.vec = sites.alongfams = NULL
        if (type=="alleles") 
        {
            minor.allele.alongfams = NULL
            for (i in 1:length(sites))
            {
                fams.site = unique(ped.mat[ped.mat[,6]==2 & (ped.mat[,5+2*sites[i]]==minor.allele.vec[i] | ped.mat[,6+2*sites[i]]==minor.allele.vec[i]),1])
                if (is.factor(fams.site)) fams.site=as.character(fams.site)
                fams.vec = c(fams.vec,fams.site)
                sites.alongfams = c(sites.alongfams,rep(sites[i],length(fams.site)))
                minor.allele.alongfams = c(minor.allele.alongfams,rep(minor.allele.vec[i],length(fams.site)))
            }
        }
        else
        {
            for (i in 1:length(sites))
            {
                # Remove subjects with missing genotype
                ped.obs = ped.mat[!is.na(ped.mat[,6+sites[i]]),]
                fams.site = unique(ped.obs[ped.obs[,6]==2 & ped.obs[,6+sites[i]]>0,1])
                if (is.factor(fams.site)) fams.site=as.character(fams.site)
                fams.vec = c(fams.vec,fams.site)
                sites.alongfams = c(sites.alongfams,rep(sites[i],length(fams.site)))            
            }
        }
    }
    else 
    {
        if (length(sites)!=length(fams)) stop ("Lengths of fams and sites vectors differs.")
        fams.vec = fams
        sites.alongfams = sites
        if (type=="alleles") minor.allele.alongfams = minor.allele.vec
    }
            
    fams.vec = as.character(fams.vec)
    
    famu = unique(fams.vec)
    famRVprob = famNcarriers = rep(NA,length(famu))
    names(famRVprob) = names(famNcarriers) = famu
    # Loop over the families
    for (f in 1:length(fams.vec))
    {
        # get carriers list
        if (type=="alleles")
            carriers = extract_carriers(ped.mat,sites.alongfams[f],fams.vec[f],type="alleles",minor.allele.alongfams[f])
        else carriers = extract_carriers(ped.mat,sites.alongfams[f],fams.vec[f],type=type)
                
        # Computation of RV sharing probability
        if (length(carriers)>0) 
        {
            if (fams.vec[f] %in% names(precomputed.prob)) 
                tmp = precomputed.prob[[fams.vec[f]]][length(carriers)]
            else tmp = suppressMessages(RVsharing(ped.listfams[[fams.vec[f]]],carriers=carriers))
            # If the RV has lower sharing probability, we keep it for this family
            if (is.na(famRVprob[fams.vec[f]]) || tmp < famRVprob[fams.vec[f]])
            {
                famRVprob[fams.vec[f]] = tmp
                famNcarriers[fams.vec[f]] = length(carriers)
            }
        }
    }
    # Identify number of informative families
    fam.info = names(famRVprob)[!is.na(famRVprob)]
    nfam.info = length(fam.info)
        
    if (nfam.info>1)
    {
        # Computing potential p-value
        potentialp = prod(pshare.vec[fam.info])
        # Extraction of the number of affected subjects in each informative family
        ped.info = ped.mat[ped.mat[,1]%in%fam.info,]
        maxN = tapply(ped.info[,6],ped.info[,1],function(vec) sum(vec==2))
        # Computing p-value
        # Families where variant is not shared by all affected subjects
        not = fam.info[famNcarriers[fam.info]<maxN]
        if (length(not)>0)
        {
            if (2^nfam.info <= maxdim)
            {
                pshare = list(ped.tocompute.vec=fam.info,pshare=pshare.vec[fam.info])
                pall = suppressWarnings(get.psubset(fam.info,not,pshare))
            }
            else pall = NA
        }
        else pall = potentialp
    }
    list(pall=pall,potentialp=potentialp)
}
