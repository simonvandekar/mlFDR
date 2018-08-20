#' Computes FDR adjusted p-values using a semiparametric bootstrap
#'
#' This function computes false discovery rate adjusted p-values using a
#' semiparametric bootstrap procedure to estimate the null distribution.
#' @param y An nXm matrix of data
#' @param statfun A function used to compute test statistics from y. Defaults to one sample t-statistic.
#' @param formred formula or character that can be coerced to a formula or a
#' design matrix for reduced model. If robust=TRUE then this must have one less
#' column than X.
#' @param mask File name for a binary mask file or niftiImage object.
#' @param data R data frame containing variables in form. If form and formred
#' are matrices then this can be NULL.
#' @param nboot Number of bootstrap sample used to estimate null distribution.
#' @param ... Arguments to pass to statfun.
#' @keywords bootstrap, FDR, false discovery rate
#' @importFrom logspline logspline dlogspline
#' @return Returns a list with the following values:
#' \describe{
#'   \item{p}{A vector of adjusted p-values that control the FDR.}
#'   \item{pi0}{The estimated proportion of true null tests.}
#' }
#' @export
bootFDR = function(y, statfun = function(y){
  out = list('mean'=colMeans(y), 'sd'=apply(y, 2, sd) )
  },
  nboot=10, ...){
  # sample size
  n = nrow(y)
  # number of tests
  m = ncol(y)
  
  # compute observed data values
  z = statfun(y)
  
  # bootstrap null estimation
  z0b = replicate(nboot, {z0 = statfun(y[sample(n, replace=TRUE),], ...);
                          z0 = (z0$mean - z$mean)/z0$sd;
                          z0} )
  z = z$mean/z$sd
  
  # small values are more "significant"
  z0b = -abs(c(z0b))
  z = -abs(c(z))
  
  # empirical CDF
  lsz = logspline(z, ubound=0)
  lsz0b = logspline(z0b, ubound=0)
  
  # estimated densities evaluated at observed data
  f0hat = dlogspline(z, lsz0b)
  fhat = dlogspline(z, lsz)
  
  # estimate proportion of true nulls, pi0
  pi0dat = pi0est(z, f0hat, fhat)
  pi0 = pi0dat[[1]][nrow(pi0dat[[1]]), 'pi0']

  # ordering for FDR
  o = order(z, decreasing=TRUE)
  ro = order(o)
  
  # p-values
  p = pmin(1, 2 * pnorm(z[o]) * pi0 / (m:1) * m )[ro]
  
  # output
  return(list(p=p, pi0=pi0))
}