gene.drop.fn <- function(g1,g2)
{
    if (is.na(g1) | is.na(g2))
    {
        return(NA)
    }

    index = g1 * 3 + g2 + 1

    if (index < 1 | index > 9) return ('error')

    switch(index,
        {return (0)},                                                   #0-0
        {return (sample(x = 0:1, size = 1, prob = rep(0.5,2)))},        #0-1
        {return (1)},                                                   #0-2
        {return (sample(x = 0:1, size = 1, prob = rep(0.5,2)))},        #1-0
        {return (sample(x = 0:2, size = 1, prob = c(0.25,0.5,0.25)))},  #1-1
        {return (sample(x = 1:2, size = 1, prob = rep(0.5,2)))},        #1-2
        {return (1)},                                                   #2-0
        {return (sample(x = 1:2, size = 1, prob = rep(0.5,2)))},        #2-1
        {return (2)},                                                   #2-2
    )
    return(goff)
}

GeneDropSim.fn <- function(trio.list, id, dt.vec, fd.indices, carriers=dt.vec, n = 1e3, k = 10, nf = 1)
{
    n.bail <- k*n; i <- 1;
    share.vec <- logical(n.bail); occur.vec <- logical(n.bail);
    geno.vec = rep(NA,length(id))
    names(geno.vec) = id
    geno.vec[fd.indices] = 0

    while( sum(occur.vec) < n & i <= n.bail )
    {
      founder <- sample(fd.indices,nf,replace = FALSE)
      geno.vec[founder] <- 1
      geno.vec.sim = geno.vec

      for (j in 1:length(trio.list))
        geno.vec.sim <- GeneDrop(trio.list[[j]], geno.vec.sim)

      if( any(is.na(geno.vec.sim))) stop("GeneDrop returned NA genotype.")

      share.vec[i] <- all( geno.vec.sim[carriers]==1 ) & all( geno.vec.sim[setdiff(dt.vec,carriers)]==0 )      
      occur.vec[i] <- any( geno.vec.sim[dt.vec]==1 )
      geno.vec[founder] <- 0
      i <- i + 1
    }
    return(sum(share.vec)/sum(occur.vec)) # note that these vectors have length n.bail
}

GeneDropSim.allsubsets.fn <- function(trio.list, id, dt.vec, fd.indices, carriers=dt.vec, n = 1e3, k = 10, nf = 1)
{
    n.bail <- k*n; i <- 1;
    share.vec <- character(n.bail); occur.vec <- logical(n.bail);
    geno.vec = rep(NA,length(id))
    names(geno.vec) = id
    geno.vec[fd.indices] = 0
    carriers = sort(carriers)

    while( sum(occur.vec) < n & i <= n.bail )
    {
        founder <- sample(fd.indices,nf,replace = FALSE)
        geno.vec[founder] <- 1
        geno.vec.sim = geno.vec
        for (j in 1:length(trio.list))
            geno.vec.sim <- GeneDrop(trio.list[[j]], geno.vec.sim)
        if( any(is.na(geno.vec.sim))) stop("GeneDrop returned NA genotype.")
        if (all( geno.vec.sim[setdiff(dt.vec,carriers)]==0 ) & any( geno.vec.sim[dt.vec]>0 )) share.vec[i] <- paste(carriers[geno.vec.sim[carriers]>0],collapse=";")     

        occur.vec[i] <- any( geno.vec.sim[dt.vec]>0 )
        geno.vec[founder] <- 0
        i <- i + 1
    }
    return(table(share.vec[share.vec!=""])/sum(occur.vec)) # note that these vectors have length n.bail
}

GeneDropSimExcessSharing.fn <- function(trio.list, id, dt.vec, fd.indices, phihat, RVfreq, carriers=dt.vec, ord=5, n = 1e3, k = 10)
{
    if (ord > 5) stop ("Order of the polynomial approximation cannot exceed 5.")
    if (ord <= 0) stop ("Order of the polynomial approximation must be at least 1.")
    if (any(phihat < 0)) stop ("Mean kinship between founders must be non-negative.")

    n.bail <- k*n; 
    geno.vec = rep(NA,length(id))
    names(geno.vec) = id
    geno.vec[fd.indices] = 0
  
    nrep = length(phihat)
    papprox = pMC = PFU.vec = rep(NA,nrep)
    nf = length(fd.indices)
  
    # Vector of phia
    phi.vec = compute.phi.vec(nf,2*nf-ord)

    # Loop over the estimates of phi in the vector phihat 
    for (r in 1:nrep)
    {
        theta = max(infer.theta(phihat[r],phi.vec))
        # Debugging code
        # print (theta)
        if (theta>=0)
        {
            i <- 1;
            share.vec <- logical(n.bail); occur.vec <- logical(n.bail);

            # Distribution of number of duplicated alleles
            distri = c(1,theta,theta^2/2,theta^3/6,theta^4/24,theta^5/120)[1:(ord+1)]
            distri = distri/sum(distri)

            while (sum(occur.vec) < n & i <= n.bail)
            {
                # Define vector of indices of the allele copies in the founders
                copyindex = 1:(2*nf)

                # Sample number of duplicated alleles
                nduplic = sample(0:ord,1,prob=distri)

                # Debugging code
                # print (nduplic)

                # If no frequency is specified for the RV, it implies a single copy of the RV
                if (missing(RVfreq))
                {
                    RVcopies = 1
                }
                else # Else the number of copies of the RV is sampled from a truncated binomial distribution from 1 to 2*nf-duplic
                {
                    distri.RVcopies = dbinom(1:(2*nf-nduplic),2*nf-nduplic,RVfreq)/(1-dbinom(0,2*nf-nduplic,RVfreq))
                    RVcopies = sample(1:(2*nf-nduplic),prob=distri.RVcopies)
                }
            }

            # Sample RV among the alleles
            RV = sample.int(2*nf-nduplic,RVcopies)

            # If RV is among the duplicated alleles, sample two gene copies for the RV
            if (RV[1] <= nduplic)
            RVi = sample.int(length(copyindex),2,replace = FALSE)
            # Else sample a single gene copy for the RV
            else RVi = sample.int(length(copyindex),1)
            founder <- rep(fd.indices,2)[copyindex[RVi]]
            copyindex = copyindex[-RVi]

            # If more than one copy of the RV is introduced
            if (length(RV) > 1)
            {
                for (k in 2:length(RV))
                {  
                    # If RV is among the duplicated alleles, sample two founders who introduce it
                    if (RV[k] <= nduplic)
                    RVi = sample.int(length(copyindex),2,replace = FALSE)
                    # Else sample a single founder introducing the RV
                    else RVi = sample.int(length(copyindex),1)
                    founder <- c(founder,rep(fd.indices,2)[copyindex[RVi]])
                    copyindex = copyindex[-RVi]        
                }
            }
            tab.founder = table(founder)
            founder = unique(founder)
            # Debugging code
            #print (tab.founder)
            geno.vec[founder] = tab.founder

            geno.vec.sim = geno.vec
            for (j in 1:length(trio.list))
                geno.vec.sim <- GeneDrop(trio.list[[j]], geno.vec.sim)
            if (any(is.na(geno.vec.sim))) stop("GeneDrop returned NA genotype.")
            share.vec[i] <- all( geno.vec.sim[carriers]==1 & all( geno.vec.sim[setdiff(dt.vec,carriers)]==0 ) )
            occur.vec[i] <- any( geno.vec.sim[dt.vec]==1 )
            geno.vec[founder] <- 0
            # Debugging code
            #print (geno.vec)
            i <- i + 1
        }
        pMC[r] = sum(share.vec)/sum(occur.vec) # note that these vectors have length n.bail
    }
    return(pMC)
}

get.fd.indices <- function(ped)
{
    dv = kindepth(ped)
    which(dv==0)
}
  
