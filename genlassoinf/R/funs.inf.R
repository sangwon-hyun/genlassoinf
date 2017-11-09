# Makeshift replacement of all usages of pval.fl1d to poly.pval;
# todo: change simulation code directly
pval.fl1d <- function(y, G, dik, sigma, approx=T, threshold=T, approxtype = c("gsell","rob"), u = rep(0,nrow(G))){
  ## return(poly.pval(y, G, u, dik, sigma, bits=NULL)$pv)
    stop("you're using pval.fl1d! Swtich to poly.pval")
}

##' Main p-value function
##' @export
poly.pval <- function(y, G, u, v, sigma, bits=NULL) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)

  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))

  pv = tnorm.surv(z,0,sd,vlo,vup,bits)
  return(list(pv=pv,vlo=vlo,vup=vup))
}

##' Main confidence interval function
##' @export
poly.int <- function(y, G, u, v, sigma, alpha, gridrange=c(-100,100),
                     gridpts=100, griddepth=2, flip=FALSE, bits=NULL) {

  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)

  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))

  xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  fun = function(x) { tnorm.surv(z,x,sd,vlo,vup,bits) }

  ## int = grid.search(xg,fun,alpha/2,1-alpha/2,gridpts,griddepth)
  int = grid.search(xg,fun,alpha,Inf,gridpts,griddepth) ## Added for one sided intervals
  tailarea = c(fun(int[1]),1-fun(int[2]))

  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }

  return(list(int=int,tailarea=tailarea))
}

##############################

# Assuming that grid is in sorted order from smallest to largest,
# and vals are monotonically increasing function values over the
# grid, returns the grid end points such that the corresponding
# vals are approximately equal to {val1, val2}

grid.search <- function(grid, fun, val1, val2, gridpts=100, griddepth=2) {
  n = length(grid)
  vals = fun(grid)

  ii = which(vals >= val1)
  jj = which(vals <= val2)
  if (length(ii)==0) return(c(grid[n],Inf))   # All vals < val1
  if (length(jj)==0) return(c(-Inf,grid[1]))  # All vals > val2
  # RJT: the above logic is correct ... but for simplicity, instead,
  # we could just return c(-Inf,Inf)

  i1 = min(ii); i2 = max(jj)
  if (i1==1) lo = -Inf
  else lo = grid.bsearch(grid[i1-1],grid[i1],fun,val1,gridpts,
         griddepth-1,below=TRUE)
  if (i2==n) hi = Inf
  else hi = grid.bsearch(grid[i2],grid[i2+1],fun,val2,gridpts,
         griddepth-1,below=FALSE)
  return(c(lo,hi))
}

# Repeated bin search to find the point x in the interval [left, right]
# that satisfies f(x) approx equal to val. If below=TRUE, then we seek
# x such that the above holds and f(x) <= val; else we seek f(x) >= val.

grid.bsearch <- function(left, right, fun, val, gridpts=100, griddepth=1, below=TRUE) {
  n = gridpts
  depth = 1

  while (depth <= griddepth) {
    grid = seq(left,right,length=n)
    vals = fun(grid)

    if (below) {
      ii = which(vals >= val)
      if (length(ii)==0) return(grid[n])   # All vals < val (shouldn't happen)
      if ((i0=min(ii))==1) return(grid[1]) # All vals > val (shouldn't happen)
      left = grid[i0-1]
      right = grid[i0]
    }

    else {
      ii = which(vals <= val)
      if (length(ii)==0) return(grid[1])   # All vals > val (shouldn't happen)
      if ((i0=max(ii))==n) return(grid[n]) # All vals < val (shouldn't happen)
      left = grid[i0]
      right = grid[i0+1]
    }

    depth = depth+1
  }

  return(ifelse(below, left, right))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector

tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)

  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1

  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits)

  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)

  p[o] = pp
  return(p)
}

##' Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, using
##' multi precision floating point calculations thanks to the Rmpfr package
mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                      (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }

  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)

  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)

  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
         (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

# Return Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# Riemann approximation tricks, by Max G'Sell

gsell.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  return(max.approx.frac(a/sd,b/sd,z/sd,mean/sd))
}


##############################

forwardStop <- function(pv, alpha=.10){
  if (alpha<0 || alpha>1) stop("alpha must be in [0,1]")
  if (min(pv,na.rm=T)<0 || max(pv,na.rm=T)>1) stop("pvalues must be in [0,1]")
  val=-(1/(1:length(pv)))*cumsum(log(1-pv))
  oo = which(val <= alpha)
  if (length(oo)==0) out=0
  else out = oo[length(oo)]
  return(out)
}

##############################

aicStop <- function(x, y, action, df, sigma, mult=2, ntimes=2) {
  n = length(y)
  k = length(action)
  aic = numeric(k)
  G = matrix(0,nrow=0,ncol=n)
  u = numeric(0)
  count = 0

  for (i in 1:k) {
    A = action[1:i]
    aic[i] = sum(lsfit(x[,A],y,intercept=F)$res^2) + mult*sigma^2*df[i]

    j = action[i]
    if (i==1) xtil = x[,j]
    else xtil = lsfit(x[,action[1:(i-1)]],x[,j],intercept=F)$res
    s = sign(sum(xtil*y))

    if (i==1 || aic[i] <= aic[i-1]) {
      G = rbind(G,s*xtil/sqrt(sum(xtil^2)))
      u = c(u,sqrt(mult)*sigma)
      count = 0
    }

    else {
      G = rbind(G,-s*xtil/sqrt(sum(xtil^2)))
      u = c(u,-sqrt(mult)*sigma)
      count = count+1
      if (count == ntimes) break
    }
  }

  if (i < k) {
    khat = i - ntimes
    aic = aic[1:i]
  }
  else khat = k

  return(list(khat=khat,G=G,u=u,aic=aic,stopped=(i<k)))
}

#these next two functions are used by the binomial and Cox options of fixedLassoInf

mypoly.pval.lee=
function(y, A, b, eta, Sigma, bits=NULL) {
    # compute pvalues from poly lemma:  full version from Lee et al for full matrix Sigma
    nn=length(y)
    eta=as.vector(eta)
  temp = sum(eta*y)
   vv=as.numeric(matrix(eta,nrow=1,ncol=nn)%*%Sigma%*%eta)
   cc = Sigma%*%eta/vv

 z=(diag(nn)-matrix(cc,ncol=1)%*%eta)%*%y
    rho=A%*%cc

  vec = (b- A %*% z)/rho
  vlo = suppressWarnings(max(vec[rho<0]))
  vup = suppressWarnings(min(vec[rho>0]))
  sd=sqrt(vv)
  pv = tnorm.surv(temp,0,sd,vlo,vup,bits)
  return(list(pv=pv,vlo=vlo,vup=vup,sd=sd))
}



mypoly.int.lee=
   function(y,eta,vlo,vup,sd, alpha, gridrange=c(-100,100),gridpts=100, griddepth=2, flip=FALSE, bits=NULL) {
    # compute sel intervals from poly lemmma, full version from Lee et al for full matrix Sigma

  temp = sum(eta*y)

  xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  fun = function(x) { tnorm.surv(temp,x,sd,vlo,vup,bits) }

  int = grid.search(xg,fun,alpha/2,1-alpha/2,gridpts,griddepth)
  tailarea = c(fun(int[1]),1-fun(int[2]))

  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }

  return(list(int=int,tailarea=tailarea))
}



mydiag=function(x){
    if(length(x)==1) out=x
    if(length(x)>1) out=diag(x)
       return(out)
   }



################################################
## Added from binSegInf R package, in https://github.com/linnylin92/binSegInf,
## as of Sep 6th, because the package is not publicly available (yet)
#' TG confidence interval from post-selection inference. Assumes that one-sided
#' test is in the right-sided direction i.e. H_0: v^T\mu>0.
#'
#' @param y numeric vector
#' @param polyhedra polyhedra object, or just a list that contains a numeric
#'     matrix |gamma| and vector |u|.
#' @param contrast contrast numeric vector
#' @param sigma numeric to denote the sd of the residuals
#' @param alpha numeric between 0 and 1 with default of 0.95. This is the
#'     significance level, guaranteeing that alpha percentage of the intervals
#'     will cover the true parameter.
#' @param gridsize numeric to denote how fine of a grid to invert the hypothesis
#'     test
#' @param alternative string of either "one.sided" or "two.sided" for the
#'     alternative. If one.sided, the alternative means the test statistic is
#'     positive.
#' @param precBits precision of Rmpfri
#' @param fac A positive factor that makes the search over \eqn{v^T\mu} to be
#'     over a range of \code{(-fac*(max(y)-min(y)), fac*(max(y)- min(y)))}.
#'     Defaults to \code{10}.
#'
#' @return a vector of two numbers, the lower and upper end of the confidence interval
#' @export
confidence_interval <- function(y, polyhedra, contrast, sigma = 1, alpha = 0.05,
                                gridsize = 250, alternative = c("two.sided", "one.sided"), precBits = NA,
                                fac=10){
    alternative <- match.arg(alternative, c("two.sided", "one.sided"))
    coverage = 1-alpha

    diff <- max(y) - min(y)
    seq.val <- seq(-fac*diff, fac*diff, length.out = gridsize)
    myfun <- get_nonzero_mean_pv_fun(y,polyhedra$gamma, polyhedra$u, contrast,sigma,alpha)
    pvalue <- sapply(seq.val, myfun)

    if(alternative == "two.sided"){

        idx <- c(.select_index(pvalue, alpha/2, T),
                 .select_index(pvalue, 1-alpha/2, F))
        c(seq.val[idx[1]], seq.val[idx[2]])

    } else if (alternative == "one.sided"){

        idx <- .select_index(pvalue, alpha, T)
        c(seq.val[idx], Inf)

    } else {
        stop("|alternative| option not known.")
    }
}


## TODO: when committing, this goes in as a commit message:
## Just so you know (Kevin), I'm replacing parts of your confidence interval function, especially the part that inserts get_nonzero_mean_pv_fun()
## - There is about a 4-5 times speedup.
## - You wrote an option called |alpha| that actually means 1-\alpha, if we follow the convention that \alpha to be the significance level of a test.
## - The grid of $v^T\mu$ you were searching for was too small; I added an option called |fac|, which defaults to 10 that allows for c(-fac*diff, fac*diff), instead of the c(-2*diff, 2*diff) that you originally wrote.



##' Selects the index of \code{vec} that is larger than \code{alpha}. Or
##' something like that.
.select_index <- function(vec, alpha, lower = T){
    idx <- ifelse(lower, min(which(vec >= alpha)), max(which(vec <= alpha)))
    if(length(idx) == 0 | is.na(idx) | is.infinite(idx)){
        warning("numeric precision suspected to be too low")
        if(lower) return(1) else return(length(vec))
    }

    if(lower & vec[idx] > alpha & idx > 1) idx <- idx - 1
    if(!lower & vec[idx] < alpha & idx < length(vec)) idx <- idx + 1

    idx
}



##' Produces a function whose single input is the null value of \eqn{v^T\mu},
##' and whose output is the TG p-value according to that nonzero null.
get_nonzero_mean_pv_fun <- function(y, G, u, v, sigma, alpha,
                                    gridrange=c(-100,100),gridpts=100, bits=NULL) {

    z = sum(v*y)
    vv = sum(v^2)
    sd = sigma*sqrt(vv)

    rho = G %*% v / vv
    vec = (u - G %*% y + rho*z) / rho
    vlo = suppressWarnings(max(vec[rho>0]))
    vup = suppressWarnings(min(vec[rho<0]))

    xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
    fun = function(x) { tnorm.surv(z,x,sd,vlo,vup,bits) }

    return(fun)
}



##' A more barebones, readable version of a one-sided CI from a
##' right-sided-alternative TG test. Matches the output of confidence_interval()
##' when using the \code{alternative="one-sided"} option.
my.one.sided.ci <- function(y,poly,contrast,sigma,alpha,fac=10,
                            gridsize=250, griddepth=2){
    ## Basic check
    if(sum(contrast*y)<0) warning("Is your contrast designed for a right-direction one-sided
test? I doubt it, since you currently have a negative sum(contrast*y)..")

    ## Make a range of v^T\mu to scan over.
    diff <- max(y) - min(y)
    seq.val <- seq(-fac*diff, fac*diff, length.out = gridsize)

    ## For each value of hypothetical values of v^TY, calculate one-sided test
    ## result
    pvs = binSegInf::pvalue(y, poly, contrast, sigma = sigma, null_mean = seq.val,
                 alternative = c("one.sided"), precBits = NA)

    ## Reject
    reject = (pvs < alpha)
    cutoff = seq.val[min(which(!reject))]

    ## Calculate the one sided confidence interval
    return(c(cutoff, +Inf))
}


##' Makes a |polyhedra| class object, which is just a list that contains two
##' elements names |gamma| and |u|, which are the matrix and vector
##' characterizing the polyhedron
##' @export
polyhedra_from_genlasso <- function(obj, numSteps){
        return(structure(list(gamma=obj$Gobj.naive$G[1:obj$nkstep[numSteps],],
                              u = obj$Gobj.naive$u[1:obj$
                                                   nkstep[numSteps]]),
               class = "polyhedra"))
}
