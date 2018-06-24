## ##' Wrapper for running 1d fused lasso. Copied over from \code{binseginf}
## ##' package on Jun 8th, 2018. 
## ##' @param invalid temporary addition
## ##' @export
## fl1d <- function(y, numSteps, sigma.add=NULL, ic.stop=FALSE, maxSteps=NULL, invalid=FALSE){

##     ## Basic checks
##     if(ic.stop){
##         ## warning("|numSteps| option will be ignored and replaced with maxSteps.")
##         assert_that(!is.null(maxSteps))
##         numSteps = maxSteps
##     }    
##     if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")
##     if(numSteps <= 0) step("numSteps must be at least 1.")
    

##     if(!is.null(sigma.add)){
##         y.addnoise = rnorm(length(y), 0, sigma.add)
##         y = y + y.addnoise
##     }

##     obj = dualpathSvd2(y, maxsteps=numSteps,
##                        D=makeDmat(length(y), ord=0), invalid=invalid)
##     obj$numSteps = numSteps

##     ## If applicable, collect IC stoppage
##     obj$ic.stop = ic.stop
##     if(ic.stop){

##         ## Obtain IC information
##         ic_obj = get_ic(obj$cp, obj$y, 2, sigma)
##         obj$stoptime = ic_obj$stoptime
##         obj$consec = consec
##         obj$ic_poly = ic_obj$poly
##         obj$ic_flag = ic_obj$flag
##         obj$ic_obj = ic_obj

##         ## Update changepoints with stopped model
##         obj$cp.all = obj$cp
##         obj$cp.sign.all = obj$cp.sign
##         if(ic_obj$flag=="normal"){
##             obj$cp = obj$cp[1:obj$stoptime]
##             obj$cp.sign = obj$cp.sign[1:obj$stoptime]
##         } else {
##             obj$cp = obj$cp.sign = c()
##         }
##     }

##     if(!is.null(sigma.add)){
##         obj$noisy = TRUE
##         obj$sigma.add = sigma.add
##         obj$y.addnoise = y.addnoise
##     }
##     return(obj)
## }

## fl1d <- function(y, numSteps){
##     D = makeDmat(m=length(y), ord=0)
##     genlassoinf(y, D=D, maxsteps=numSteps)
## }

##' Generic interface function. Partially ported over from |genlasso| package.
##' @param type One of \code{c("fl1d", "fl2d", "linear.tf")}, which calls some
##'     convenience functions. Defaults to NULL.
genlassoinf <- function(y, X, D=NULL, type=NULL, approx=FALSE, maxsteps=2000, minlam=0,
                        rtol=1e-7, btol=1e-7, eps=1e-4, verbose=FALSE,
                        svd=TRUE) {
 
    ## Basic checks
    if (missing(y)) stop("y is missing.")
    if (!is.numeric(y)) stop("y must be numeric.")
    if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 1].")
    if (missing(X)) X = NULL
    if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
    if(is.null(D) ){
        if(is.null(type)){
            stop("Must provide either D or type")
        } else {
            if(type == "fl1d"){
                D = makeDmat(length(y), ord=0)
            } else if (type == "linear.tf"){
                D = makeDmat(length(y), ord=1)
            } else if (type == "fl2d"){
                D = makeDmat(length(y), type="2d")
            } else {
                stop("specialty |type| does not exist yet!")
            }
        }
    }
    if (!is.null(D)){
        if (!is.matrix(D) && c(attributes(class(D))$package,"")[[1]] != "Matrix") {
          stop("D must be a matrix or a Matrix (from the Matrix package).")
        }
        if (is.null(X) && length(y)!=ncol(D)) stop("Dimensions don't match [length(y) != ncol(D)].")
        if (checkrows(D)) stop("D cannot have duplicate rows.")
    }


  # For simplicity
  y = as.numeric(y)
  
  # X should be treated as the identity
  if (is.null(X)) {
      out = dualpathSvd(y,D,approx,maxsteps,minlam,rtol,btol,verbose)
  } else {
      stop("Not written yet! Although the treatment from genlasso::genlasso() is already written in code and should work.")
      if (!is.matrix(D)) {
          warning("Converting D to a dense matrix, because X is not the identity.")
          D = as.matrix(D)
      }
         n = nrow(X)
      p = ncol(X)
      if (length(y)!=n) stop("Dimensions don't match [length(y) != nrow(X)].")
      if (ncol(D)!=p) stop("Dimensions don't match [ncol(X) != ncol(D)].")
         ridge = FALSE
      if (p > n) {
          if (eps<=0) stop("eps must be positive when X has more columns than rows.")
          warning(sprintf("Adding a small ridge penalty (multiplier %g), because X has more columns than rows.",eps))
          ridge = TRUE
      } else {
          ## Check that X has full column rank
          x = svd(X)
          if (all(x$d >= rtol)) {
              y2 = as.numeric(x$u %*% t(x$u) %*% y)
              Xi = x$v %*% (t(x$u) / x$d)
              D2 = D %*% Xi
          }
          else {
              if (eps<=0) stop("eps must be positive when X is column rank deficient.")
              warning(sprintf("Adding a small ridge penalty (multiplier %g), because X is column rank deficient.",eps))
              ridge = TRUE
          }
      }
      
      if (ridge) {
          x = svd(rbind(X,diag(sqrt(eps),p)))
          y2 = as.numeric(x$u %*% t(x$u) %*% c(y,rep(0,p)))
          Xi = x$v %*% (t(x$u) / x$d)
          D2 = D %*% Xi
      }
      
      if (!svd) stop("Only svd option is written so far!") ##out = dualpath(y2,D2,approx,maxsteps,minlam,rtol,btol,verbose)
      else out = dualpathSvd2(y2,D2,approx,maxsteps,minlam,rtol,btol,verbose)
      
      ## Save these path objects for internal use later
      out$pathobjs$y2 = y2
      out$pathobjs$Xi = Xi
         # Fix beta, fit, y, bls, and save the X matrix
      out$beta = Xi %*% out$fit
      out$fit = X %*% out$beta
      out$y = y
      out$bls = Xi %*% y2
      out$X = X
  
      # Fix df
      if (ridge) out$df = out$df - n
      else out$df = out$df - (n-p)
  }
    out$call = match.call()
    class(out) = c("genlassoinf", "list")
    return(out)
}
