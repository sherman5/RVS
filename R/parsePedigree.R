# parses all pedigree information into a simple list of Subject objects
parsePedigree <- function(ped)
{
    # create empty list of subjects
    subjects <- replicate(length(ped[['id']]), new('Subject'))

    # get info for each subject
    for (i in 1:length(subjects))
    {
        # get mom and dad index
        subjects[[i]]@mom <- ped[['mindex']][[i]]
        subjects[[i]]@dad <- ped[['findex']][[i]]

        # find all subjects that refer to this one as 'mom' or 'dad'    
        subjects[[i]]@offspring <- c(which(ped[['findex']] == i),
            which(ped[['mindex']] == i))

        # find spouses
        for (c in subjects[[i]]@offspring)
        {
            subjects[[i]]@spouses <- c(subjects[[i]]@spouses,
                ped[['mindex']][[c]], ped[['findex']][[c]])
        }

        # remove duplicates and self
        subjects[[i]]@spouses <- unique(subjects[[i]]@spouses)
        subjects[[i]]@spouses <- 
            subjects[[i]]@spouses[subjects[[i]]@spouses != i]        
    }
    return (subjects)
}

# get the indices of the founders
getFounders <- function(subjects)
{
    sumParents <- function(sub) {return (sub@mom + sub@dad)}
    sumInLaws <- function(sub)
        {return (sum(unlist(sapply(subjects[sub@spouses], sumParents))))}

    parents <- sapply(subjects, sumParents)
    inLaws <- sapply(subjects, sumInLaws)

    return (which(parents == 0 & inLaws == 0))
}

# get the indices of the final descendants
getFinalDescendants <- function(subjects)
{
    num <- sapply(subjects, function(sub) {return (length(sub@offspring))})
    return (which(num == 0))
}

# get indices of all children
getNextGeneration <- function(subjects)
{
    allChildren <- sapply(subjects, function(sub) {return (sub@offspring)}) 
    return (unique(unlist(allChildren)))
}

# get a list of all descendants from each founder
getDescendants <- function(founders, subjects, ped)
{
    descendants <- list()
    for (f in 1:length(founders))
    {
        # get first generation after this founder
        gen <- getNextGeneration(subjects[founders[f]])
        descendants[[f]] <- sapply(gen, function(x) {return (c(x,1))})

        # check if no descendants
        if (length(gen) == 0)
        {
            stop("error, founder ", ped[['id']][founders[f]],
                " has no descendants")
        }

        # get all subsequent generations
        depth <- 2
        gen <- getNextGeneration(subjects[gen])
        while (length(gen) > 0)
        {
            descendants[[f]] <- cbind(descendants[[f]], sapply(gen,
                function(x) {return (c(x,depth))}))
            
            depth <- depth + 1
            gen <- getNextGeneration(subjects[gen])
        }
    }
    return (descendants)
}
        
# compute RV sharing probability
RVsharingNew <- function(ped)
{
    subjects <- parsePedigree(ped)
    founders <- getFounders(subjects)
    descendants <- getDescendants(founders, subjects, ped)

    numer = 0
    denom = 0

    for (f in 1:length(founders))
    {
        dist <- descendants[[f]][2, descendants[[f]][1,] > 0]

        commonAncestor <- 1

        numer <- numer + 0.5 ^ sum(dist) * commonAncestor
        denom <- denom + 1 - prod(1 - 0.5 ^ dist)
    }
    return (numer / denom)
}


