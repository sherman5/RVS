listRVhaplo = function(vec,chr,legend.dat,haplo.dat,filteredvcf.dat,sample.dat,selected.dat,CV.bool,build="hg19",minCV=2)
{
	# vec : line of a data frame containing [1] the gene name, [4] the line for the first and 
	#		[5] last variant in the gene in the data frames legend.dat and haplo.dat.
	# legend.dat : data frame containing an id variable with the variant identifiers and a position variable with its position.
	# haplo.dat : data frame with one column for each haplotypes of the study subjects
	# filteredvcf.dat : optional data frame containing an ID variable with the variants to be included
	# 					in the analysis and a RefgeneGeneName column with the gene name.
	# sample.dat : data.frame with a variable ID_1 identifying the family.
	# selected.dat : optional logical vector of weather each variant in legend.dat meets the desired selection criterion
	famid = as.character(sample.dat$ID_1)
         # gene name
        g = as.character(unlist(vec[1]))
        # Read section of haplotype file for the gene
        # Subtract 1 because there is no header line.
        phasedvariants = legend.dat[(as.numeric(vec[4]):as.numeric(vec[5]))-1,]
        haplotypes = haplo.dat[(as.numeric(vec[4]):as.numeric(vec[5]))-1,]
        
        if (!is.matrix(haplotypes)) haplotypes = matrix(haplotypes,1,length(haplotypes))

	  if (missing(selected.dat)) selected = rep(TRUE,nrow(phasedvariants))
	  else selected = selected.dat[(as.numeric(vec[4]):as.numeric(vec[5]))-1]
	  
	  if (missing(filteredvcf.dat)) selected2 = rep(TRUE,nrow(phasedvariants))
	  else {
	    # Extract name (for dbSNP variants) and position of the gene's RVs from vcf file
        selectedvars = as.character(filteredvcf.dat$ID[filteredvcf.dat$RefgeneGeneName==g])
        temp = strsplit(substr(selectedvars,nchar(chr)+2,nchar(selectedvars)),"_")
        selectedpos = ifelse(substr(selectedvars,1,2)=="rs",0,sapply(temp,function(v) as.numeric(v[1])))
        selected2 = phasedvariants$id%in%selectedvars | phasedvariants$position%in%selectedpos
        }
        
        # Extract RV haplotypes from haplotype file
if (any(selected&selected2))
{
haplotypes.tmp = haplotypes[selected&selected2,]

        if (!is.matrix(haplotypes.tmp)) haplotypes.tmp = matrix(haplotypes.tmp,1,length(haplotypes.tmp))

# Switch 0 and 1 when allele 1 is more frequent
RVhaplotypes = t(apply(haplotypes.tmp,1,function(vec) if (mean(vec)>0.5) 1-vec else vec))

# Variants not selected with rs number are common variants
if (missing(CV.bool)) CV.bool = substr(phasedvariants$id,1,2)=="rs" & !(selected&selected2)
if (sum(CV.bool)>=minCV) common_var_haplo = haplotypes[CV.bool ,]
# If all variants are selected, take dbSNP variants to define haplotypes
else common_var_haplo = haplotypes[substr(phasedvariants$id,1,2)=="rs" ,]

if (!is.matrix(common_var_haplo )) common_var_haplo  = matrix(common_var_haplo,1,ncol(haplotypes))
if (!is.matrix(RVhaplotypes )) RVhaplotypes  = matrix(RVhaplotypes,1,ncol(haplotypes))

# If there is at least one dbSNP variant in the gene discard RVs
# on multiple haplotypes defined by these RVs
if (nrow(common_var_haplo)>0)
{
trueRV = apply(RVhaplotypes,1,ncvhaplo,chm = common_var_haplo,fams=rep(famid,rep(2,length(famid))))==1
if (sum(trueRV)>0)
{
trueRVhaplotypes = RVhaplotypes[trueRV,]
if (!is.matrix(trueRVhaplotypes )) trueRVhaplotypes  = matrix(trueRVhaplotypes,1,ncol(haplotypes))

# Add the family ID as an additional "variant", so the same common haplotype in different families 
# is assigned a different common_var_string
tmp = rbind(common_var_haplo,rep(famid,rep(2,length(famid))))
common_var_string = apply(tmp,2,function(vec) paste(vec,collapse=""))
recodedRV.vec = apply(trueRVhaplotypes,2,sum)

# Recoding haplotypes into variants

# for each common haplotype, save the vector of number of RVs on each copy of that haplotype 
# (or NULL if no RV on the haplotype)
recodedRV.list = tapply(recodedRV.vec,common_var_string,function(vec) if(any(vec>0)) vec else NULL )
RVonhaplo = !sapply(recodedRV.list,is.null)
recodedRV.mat = matrix(0,sum(RVonhaplo),ncol(haplotypes))
# The row names are the common variant haplotypes with RVs on them (without specification of which RV)
dimnames(recodedRV.mat) = list(names(recodedRV.list)[RVonhaplo],common_var_string)
for (i in names(recodedRV.list)[RVonhaplo])
	# Each common haplotype copy present in the haplotype sample is assigned the number of RVs it carries on 
	# the recoded haplotypes (should be the same for all)
	recodedRV.mat[i,dimnames(recodedRV.mat)[[2]]==i] = recodedRV.list[[i]]

tt = table(trueRV)
cat(g,tt,"\n")
}
else recodedRV.mat = NULL
}
else 
	{
	recodedRV.mat = apply(RVhaplotypes,2,sum)
	if (!is.matrix(recodedRV.mat)) recodedRV.mat = matrix(recodedRV.mat,1,ncol(haplotypes))	
	}
}
else recodedRV.mat = NULL

recodedRV.mat
}
