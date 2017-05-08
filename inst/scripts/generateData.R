library(gRbase)

sampleDAG <- dag(1,2,3,4,5,6,7,8,9,10,11,12,13,
    c(4,1), c(4,2),
    c(7,2), c(7,3),
    c(8,4), c(8,5),
    c(9,4), c(9,5),
    c(13,6), c(13,7),
    c(10,6), c(10,7),
    c(11,8), c(11,13),
    c(12,9), c(12,10))

save(sampleDAG, file='sampleGraphs.RData')
