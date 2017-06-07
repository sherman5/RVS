#' \code{inferNumAlleles} Returns the most likely number of distinct alleles
#'  among nf founders based on mean estimated kinship phi
#' @keywords internal
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

#' \code{computePhiVec} the vector of expected phi_a for nf founders for
#'  numbers of distinct alleles a from amin to 2*nf-1
#' @keywords internal
computePhiVec <- function(nf, amin=2*nf-2)
{
    a = amin:(2*nf-1)
    term1 = (2*nf-a)*(2*nf-a-1) / (nf*(nf-1))
    term2 = (a-nf)*(2*nf-a) / (nf*(nf-1)) + 2*(a-nf)*(2*nf-a) / (nf*(2*nf-1))
    return((0.5*term1 + 0.25*term2) / (nf-1))
}

#' \code{inferTheta} solve the parameter theta for polynomial approximation
#'  of the distribution of the number of distinct alleles. This is a general
#' function for polynomials of order 2 to 5.
#'
#' @param phi the mean estimated kinship between founders
#' @param phiVec contains phi_a for a = 2*nf-ord to 2*nf-1, where
#'  ord must be between 2 and 5
#' @return real roots of the polynomial approximation
#' @keywords internal
inferTheta <- function(phi, phiVec)
{
    ord = length(phi.vec)
    phi.diff = (phi - phi.vec)
    coef.vec = 1/factorial(1:ord)
    racines = polyroot(c(phi,phi.diff[ord:1]*coef.vec))
    return(Re(racines)[abs(Im(racines))<1e-10]) # Return only the real roots
}

#' \code{computePFU} computation of P[FjU] using equation 21 of Bureay et al.
#'
#' @param nf number of founders of the pedigree
#' @param theta value of the parameter of the polynomial distribution
#' @param ord order of the polynomial approximation to the distribtion of the
#'  number of distinct alleles in the founders (noted d in Bureay et al.).
#'  Must be <= 5
#' @return P[FjU] (scalar)
#' @keywords internal
computePFU <- function(nf, theta, ord=2)
{
    a = (2*nf):(2*nf-ord)
    distri = c(1,theta,theta^2/2,theta^3/6,theta^4/24,theta^5/120)[1:(ord+1)]
    return(weighted.mean(2/nf - 2/a,distri))
}

