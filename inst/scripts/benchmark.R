library(RVsharing)
source('~/Work/RVsharing/ClonedRepo/inst/scripts/generateData.R')

# time relative to number of nodes

nGen <- 10
bm <- matrix(nrow = nGen, ncol = 2)
for (n in 1:nGen)
{
  ped <- simPedigree(n, function() 2, 0)
  len <- length(ped$affected)
  ped$affected[(len-2):(len-1)] <- 1
  bm[n,1] <- length(ped$id)  

  print(paste('running', n))
  bm[n,2] <- system.time(RVsharing(ped, alleleFreq = 1e-9))[1]
}

const <- bm[nGen,2] / (bm[nGen,1] ^ 2)
plot(bm[,1], bm[,2])
fit <- const * (bm[,1] ^ 2)
lines(bm[,1], fit)

# time relative to size of cliques

ped <- list()
ped$id <- 1:14
ped$findex <- c(0,0,rep(1,14),3,5,7,9)
ped$mindex <- c(0,0,rep(2,14),4,6,8,10)
ped$sex <- c(rbind(rep(1,7), rep(2,7)))
ped$affected <- rep(0, 14)
ped$affected[c(13,14)] <- 1
class(ped) <- 'pedigree'
plot(ped)
RVsharing(ped)

# time relative to num affected

ped <- simPedigree(6, function() 2, 0)
parents <- sapply(ped$id, function(i) c(ped$findex[i], ped$mindex[i]))
nonFounders <- which(parents[1,] != 0)

bm <- matrix(nrow=20, ncol=2)
for (nAff in seq(1,100,5))
{
  aff <- sample(nonFounders, nAff+1)
  ped$affected[aff] <- 1
  n <- ceiling(nAff / 5)
  bm[n,1] <- nAff
  bm[n,2] <- system.time(RVsharing(ped, alleleFreq = 1e-9))[1]
  ped$affected[aff] <- 0
}

plot(bm[,1], bm[,2])
