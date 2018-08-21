#' Computes pi0 via EM algorithm
#'
#' This function computes the proportion of true null hypothesis conditional on
#' an estimate of the marginal and null distributions.
#' @param z A vector of t-statistics.
#' @param f0hat A vector of length length(z) that is the null PDF evaluated at
#' z.
#' @param fhat A vector of length length(z) that is the null PDF evaluated at z.
#' @return Returns a data frame with the estimation results.
pi0est = function(z, f0hat, fhat, pi0=0.5, eps=0.001){
  maxit = 50; it = 0
  out = data.frame(iter=0:maxit, pi0=rep(NA, maxit+1), conv=rep(NA,maxit+1))
  out$pi0[1] <- pi0
  pi0t = pi0 + 2 * eps
  while(abs(pi0 - pi0t)>eps & it<=maxit ){
    pi0t = pi0
    # probabilities for the latent variables
    p_xj = pmin(1, pi0t * f0hat /fhat)
    
    #Q = function(pi0){
    # (1-p_xj) %*% log( (fhat - pi0 * f0hat) / (1-pi0) )
    #}
    #vQ = Vectorize(Q, vectorize.args = 'pi0')
    dQ = function(pi0){
      (1-p_xj) %*% (f0hat / (fhat - pi0 * f0hat) - (1-pi0)^-1 )
    }
    # upper boundary for pi0
    pi0 = pmin(1, min(fhat/f0hat))
    # estimate of unknown parameter
    root = tryCatch(uniroot(dQ, c(0,pi0)), error=function(x){ return(list(root=pi0))})
    it = it + 1
    out[ out$iter==it, 'pi0'] = root$root 
    out[ out$iter==it, 'pi0max'] = pi0
    pi0=root$root
  }
  out = out[!is.na(out$pi0),]
  out = list(out, convergence=(abs(pi0 - pi0t)<=eps) )
}
