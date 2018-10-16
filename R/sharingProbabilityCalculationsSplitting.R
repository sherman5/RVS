#' calculate sharing probability in basic case
#' @keywords internal
#'
#' @description Assume that only one founder can introduce the variant to 
#' the pedigree. Condition on each founder and sum over all resulting
#' probabilities. 
#' @param procPed pedigree that has been through processPedigree()
#' @return sharing probability
oneFounderSharingProbSplitting <- function(procPed)
{
	nf = length(procPed$founders)
	# Creating genealogy with procPed
	gen = GENLIB::gen.genealogy(cbind(procPed$id,procPed$parents,procPed$ped$sex))

	# Find occurence of founders in ancestry of each affected subject
	affByFounder = GENLIB::gen.occ(gen,pro=procPed$affected,procPed$founders)
	
    # sum over probs, conditioning on each founder introducing variant
	carrier.sets = list()
	for (i in length(procPed$affected):1)
	carrier.sets = c(carrier.sets, combn(procPed$affected,i,simplify=FALSE))
    carrier.numer <- carrier.noRV <- rep(0,length(carrier.sets))
    for (f in procPed$founders) #TODO: use sapply here
    {
    	# Extract subpedigree 
    	subaffected = procPed$affected[affByFounder[,as.character(f)]>0]
    	subgen = gen.branching(gen,pro=subaffected,ances=f)
    	subped = gen.genout(subgen)
    	# Adding dummy parents
		mb = subped$father>0 & subped$mother==0
		pb = subped$father==0 & subped$mother>0
		for (p in subped$ind[pb]) subped$father[subped$ind==p] = gen.parent(gen.famatteints,p,output="Fa")
		for (m in subped$ind[mb]) subped$mother[subped$ind==m] = gen.parent(gen.famatteints,m,output="Mo")
		p.vec = gen.parent(gen.famatteints,subped$ind[pb],output="Fa")
		m.vec = gen.parent(gen.famatteints,subped$ind[mb],output="Mo")
		pedtmp=data.frame(ind=c(p.vec,m.vec),father=0,mother=0,sex=c(rep(1,length(p.vec)),rep(2,length(m.vec))))
		subped = rbind(subped,pedtmp)
		
		# Recreating a processed ped
		subprocPed = list('parents'=cbind(subped$father,subped$mother),'id'=subped$ind,'affected'=subaffected,'founders'=which(subped$father == 0))
		
		    # set all founders to 0 (no variant)	 except current founder f	    
	    net <- try(createNetwork(subprocPed))
	    if (class(net)!="try-error")
	    {
	    evid = rep('0', length(subprocPed$founders))
	    evid[f] = '1';
	    net <- gRain::setEvidence(net, as.character(subprocPed$founders),evid)		
		
        # compute probabilities
        carrier.noRV <- carrier.noRV + sapply(carrier.sets,noRVProbVec,net=net, procPed=subprocPed)
        carrier.numer <- carrier.numer + sapply(carrier.sets,numerProbVec,net=net, procPed=subprocPed)
        }
    }
    return()
}
