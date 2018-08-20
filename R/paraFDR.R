#' Computes FDR adjusted p-values parametrically
#'
#' This function computes false discovery rate adjusted p-values using a
#' parametric null distribution (dnorm(0,1)).
#' @param z A vector of t-statistics.
#' @importFrom logspline logspline dlogspline
#' @return Returns a list with the following values:
#' \describe{
#'   \item{p}{A vector of adjusted p-values that control the FDR.}
#'   \item{pi0}{The estimated proportion of true null tests.}
#' }
#' @export
paraFDR = function(z){
  m = length(z)
  z = -abs(z)
  
  # empirical CDF
  lsz = logspline(c(z), ubound=0)
  
  # estimated densities evaluated at observed data
  f0hat = 2*dnorm(z)
  fhat = dlogspline(z, lsz)
  # compute densities
  pi0dat = pi0est(z, f0hat, fhat)
  pi0 = pi0dat[[1]][nrow(pi0dat[[1]]), 'pi0']
  
  o = order(z, decreasing=TRUE)
  ro = order(o)
  
  p = pmin(1, 2 * pnorm(z[o]) * pi0 / (m:1) * m )[ro]
  return(list(p=p, pi0=pi0))
}