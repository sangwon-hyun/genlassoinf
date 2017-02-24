##' Function to create a confidence band of a linear contrast of the mean
##' vector, from the simultaneous confidence band of the mean based on the SMUCE
##' algorithm of Munk et al. (2015)
##' @param y Data vector
##' @param d a vector of the same length of \code{y}
##' @example example/smuce-example.R
CI_smuce = function(y, d, alpha=.05){

    ## Basic checks
    stopifnot(length(y)==length(d))

    ## Form pointwise min/max, then sum
    v = d/sqrt(sum(d^2))
    sm = smuceR(y,jumpint=T,confband=T,alpha=alpha)
    cb = confband(sm)
    cb.unordered.pointwise = cb[,2:3] * v
    cb.ordered.pointwise = apply(cb.unordered.pointwise,1,function(cc)range(cc) )
    return(apply(cb.ordered.pointwise,1,sum))
}
