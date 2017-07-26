#' most likely number of distinct alleles among founders
#' @keywords internal
#'
#' @description Calculates the most likely number of distinct alleles
#'  among nf founders based on the mean estimated kinship coefficient
#' @param phi mean estimated kinship coefficient
#' @param nf number of founders
#' @return number of distinct alleles
inferNumAlleles <- function(phi, nf)
{
    a = nf:(2*nf)
    term1 = (2*nf-a)*(2*nf-a-1) / (nf*(nf-1))
    term2 = 2*(a-nf)*2*(2*nf-a) /
        (2*(a-nf)*2*(2*nf-a) + (a-nf)*(2*(a-nf)-1) + 2*(2*nf-a)*(2*nf-a-1))
    phi.vec = 0.5*term1 + 0.25*term2
    phi.diff = (phi.vec - phi)
    return(a[phi.diff>0 & c(phi.diff[-1],0)<0])
}

#' expected kinship coefficient for different number of alleles
#' @keywords internal
#'
#' @param nf number of founders
#' @param amin minimum number of distinct alleles
#' @return vector of expected phi_a for nf founders for
#'  numbers of distinct alleles from amin to 2*nf-1
computePhiVec <- function(nf, amin=2*nf-2)
{
    a = amin:(2*nf-1)
    term1 = (2*nf-a)*(2*nf-a-1)/(nf*(nf-1))
    term2 = (a-nf)*(2*nf-a)/(nf*(nf-1)) + 2*(a-nf)*(2*nf-a)/(nf*(2*nf-1))
    return((0.5*term1 + 0.25*term2)/(nf-1))
}

#' solve the parameter theta for polynomial approximation of the
#'  distribution of the number of distinct alleles.
#' @keywords internal
#'
#' @param phi the mean estimated kinship between founders
#' @param phiVec contains phi_a for a = 2*nf-ord to 2*nf-1, where
#'  ord must be between 2 and 5
#' @return real roots of the polynomial approximation
inferTheta <- function(phi, phiVec)
{
    ord = length(phiVec)
    phi.diff = (phi - phiVec)
    coef.vec = 1/factorial(1:ord)
    racines = polyroot(c(phi,phi.diff[ord:1]*coef.vec))
    return(max(Re(racines)[abs(Im(racines))<1e-10])) # Return only the maximal real roots
}

#' computation of P[FjU] using equation 21 of Bureau et al.
#' @keywords internal
#'
#' @param nf number of founders of the pedigree
#' @param theta value of the parameter of the polynomial distribution
#' @param ord order of the polynomial approximation to the distribtion of the
#'  number of distinct alleles in the founders (noted d in Bureay et al.).
#'  Must be <= 5
#' @return P[FjU] (scalar)
computePFU <- function(nf, theta, ord=5)
{
    a = (2*nf):(2*nf-ord)
	distri = theta^(0:ord)/factorial(0:ord)
	return(weighted.mean(2/nf - 2/a,distri))
}

#' make the neccesary correction for when founders have a non-zero
#'  kinship coefficient
#' @keywords internal
#'
#' @param nf number of founders
#' @param kinshipCoeff mean kinship coefficient among all founders
#' @param ord order of the polynomial approximation to the distribtion of the
#'  number of distinct alleles in the founders (noted d in Bureay et al.).
#'  Must be <= 5
#' @return weight used in probability calculation
#' @keywords internal
relatedFoundersCorrection <- function(nf, kinshipCoeff, ord=5)
{
    phiVec <- computePhiVec(nf, 2*nf-ord)
    theta <- inferTheta(kinshipCoeff, phiVec)
    PFU <- computePFU(nf, theta, ord)
    return(nf * PFU)
}

