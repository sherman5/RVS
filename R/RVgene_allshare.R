# ped.mat : pedigrees coded as in a ped file
# ped.listfams : list of pedigree objects, one object for each pedigree in ped.mat
# fams : vector of families carrying the variants listed in the corresponding position in sites
# sites : vector of the variant sites for each family in the fams vector 
# minor.allele.vec : vector of the minor alleles at each site in the sites vector
RVgene_allshare = function(ped.mat,ped.listfams,sites,fams,pshare.vec,type="alleles",minor.allele.vec,precomputed.prob=list(0),maxdim = 1e9)
{
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
			#cat (f,"\n")
			if (fams.vec[f] %in% names(precomputed.prob)) 
        		tmp = precomputed.prob[[fams.vec[f]]][length(carriers)]
			else tmp = RVsharing(ped.listfams[[fams.vec[f]]],carriers=carriers)@pshare
			# If the RV has lower sharing probability, we keep it for this family
			if (is.na(famRVprob[fams.vec[f]]) || tmp < famRVprob[fams.vec[f]])
			{
				famRVprob[fams.vec[f]] = tmp
				famNcarriers[fams.vec[f]] = length(carriers)
			}
		}
		#print(famRVprob)
	}

	# Identify number of informative families
	fam.info = names(famRVprob)[!is.na(famRVprob)]
    # print(fam.info)
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
    		  	pall = get.psubset(fam.info,not,pshare)
			}
    		else pall = NA
		}
		else pall = potentialp
	}
	list(pall=pall,potentialp=potentialp)
}
