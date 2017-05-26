library(RVsharing)
source('~/Work/RVsharing/ClonedRepo/inst/scripts/generateData.R')

# time relative to number of nodes

nGen <- 3
bm <- matrix(nrow = nGen, ncol = 2)
for (n in 1:nGen)
{
  ped <- simPedigree(n, function() 4, 0)
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
ped$findex <- c(0,0,rep(1,16),seq(from=3,by=2,length=8),seq(from=19,by=2,length=4),c(27,29),31)
ped$mindex <- c(0,0,rep(2,16),seq(from=4,by=2,length=8),seq(from=20,by=2,length=4),c(28,30),32)
len <- length(ped$findex)
ped$id <- 1:len
ped$sex <- c(rbind(rep(1,len/2), rep(2,len/2)),1)
ped$affected <- rep(0, len)
ped$affected[c(19,20)] <- 1
class(ped) <- 'pedigree'
plot(ped)
RVsharing(ped)

ped <- list()
ped$findex <- c(0,0,1,1,0,0,rep(3,6),rep(4,6),7,8,9,10,11,12)
ped$mindex <- c(0,0,2,2,0,0,rep(5,6),rep(6,6),13,14,15,16,17,18)
len <- length(ped$findex)
ped$id <- 1:len
ped$sex <- rep(1,len)
ped$affected <- rep(0, len)
ped$affected[c(3,4)] <- 1
class(ped) <- 'pedigree'
plot(ped)
RVsharing(ped)

ped <- list()
ped$findex <- c(rep(0,8),1,3,5,7,9,11)
ped$mindex <- c(rep(0,8),2,4,6,8,10,12)
len <- length(ped$findex)
ped$id <- 1:len
ped$sex <- rep(1,len)
ped$affected <- rep(0, len)
ped$affected[c(5,6)] <- 1
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
