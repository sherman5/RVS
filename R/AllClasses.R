#setClass("Trio", representation = list(
#    id = "character",
#    spouse = "character",
#    offspring = "list"
#))

setClass("Subject", representation = list(
    mom = 'numeric',
    dad = 'numeric',
    spouses = 'vector',
    offspring = 'vector'
))

#setClass("RVsharingProb", representation = list(
#    pshare="numeric",
#    iancestors="character",
#    desfounders="list",
#    id="character",
#    dad.id="character",
#    mom.id="character",
#    carriers="character"
#))

#setClassUnion("TrioOrChar", c("Trio","character"))


