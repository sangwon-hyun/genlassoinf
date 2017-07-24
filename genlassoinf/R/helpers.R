##' Make sure that no two rows of A are the same (this works with probability
##' one). Borrowed from https://github.com/statsmaths/genlasso
checkrows <- function(A) {
  b = rnorm(ncol(A))
  a = sort(A%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}

##' Helper function to obtain the effective design X matrix, by breaking X at
##' the fused lasso regression breakpoints.
get.augmented.X = function(X,breaks,return.groups = F,TT,J,group.inds){
    find.which.group = function(hit) which(sapply(1:J, function(jj) hit %in% group.inds[[jj]]))
    which.group.breaks = sapply(breaks, find.which.group)
    breaks.by.group = list()
    for(groupnum in 1:J){
      thisgroup.breaks = breaks[which.group.breaks==groupnum]
      thisgroup.breaks = thisgroup.breaks - TT*(groupnum-1) + (groupnum-1)
      breaks.by.group[[groupnum]] = thisgroup.breaks
    }

    # Function to break var into #|breaks| variables broken at |breaks|
    brk.var = function(var,breaks){
      augmented.breaks =c(0,sort(breaks),length(var))
      all.subvecs =  sapply(1:(length(augmented.breaks)-1),
                     function(ii){inds = (augmented.breaks[ii]+1):(augmented.breaks[ii+1])
                                  newvar = rep(0,length(var))
                                  newvar[inds] = var[inds]
                                  return(newvar)
                                 })
      return(all.subvecs)
    }

    X.augmented = do.call(cbind,
                          lapply(1:length(breaks.by.group),
                                  function(ii) brk.var(X[,ii], breaks.by.group[[ii]])))

     if(return.groups){ return(breaks.by.group) } else { return(X.augmented) }

  }

##' Helper function to take all neighboring-to-each-other clusters,
##' And declutter them by removing all but (rounded up) centroids
##' @export
declutter = function(coords, how.close = 1, sort=T, indexonly = F){#closeby.same.direction.are.disallowed=F

    ## n
    unsorted.coords=coords
    coords = sort(coords)

    ## error checking
    if(length(coords)<=1){
      if(length(coords)==0) cat('\n',"attempting to declutter", length(coords), "coordinates",'\n')
      return(coords)
    }

    ## get the clique memberships
    adjacent.diffs = abs(coords[1:(length(coords)-1)] - coords[2:length(coords)])

    cliq.num=1
    cliq.vec=rep(NA,length(coords))
    for(ii in 1:length(adjacent.diffs)){
      if(adjacent.diffs[ii] <= how.close){  ## used to be ==1
        cliq.vec[ii] = cliq.vec[ii+1] = cliq.num
      } else {
        cliq.num = cliq.num+1
      }
    }

    ## determine who will leave
    leavelist = c()
    for(cliq.num in unique(cliq.vec[!is.na(cliq.vec)])){
        members = which(cliq.vec == cliq.num)
        stay = round(mean(members))
        leave = members[members!=stay]
        leavelist = c(leavelist,leave)
    }
    if(length(leavelist)>=1){
      processed.coords = coords[-leavelist]
    } else {
      processed.coords = coords
    }

    if(indexonly){
      return(which(unsorted.coords %in% processed.coords))
    } else {
        if(sort){
        return(processed.coords)
      } else {
        return(unsorted.coords[unsorted.coords %in% processed.coords])
      }
    }
}


##' Gets k'th Vup or Vlo for user-input y, G, dik and sigma
getVs = function(y, G, dik, sigma, u = rep(0, nrow(G))){

  rho   <- G %*% dik / sum(dik^2)
  V <- (u - (G%*%y) + rho %*% (dik %*% y) ) / rho

  Vlo <- max(V[which(rho > 0)])
  Vup <- min(V[which(rho < 0)])

  Vs <- c(Vlo,Vup)
  names(Vs) <- c("Vlo","Vup")

  return(Vs)
}

# Returns a sequence of integers from a to b if a <= b,
# otherwise nothing. You have no idea how important this
# function is...
Seq = function(a,b) {
  if (a<=b) return(a:b)
  else return(integer(0))
}

# Returns the sign of x, with Sign(0)=1.
Sign <- function(x) {
  return(-1+2*(x>=0))
}


## A bunch of helpers
getDtfSparse = function(n,ord) {
  D = bandSparse(n, m=n, c(0,1), diagonals=list(rep(-1,n),rep(1,n-1)))
  D0 = D
  for (i in Seq(1,ord)) D = D0 %*% D
  return(D[Seq(1,n-ord-1),])
}

getDtf = function(n, ord) {
  return(as.matrix(fSparse(n,ord)))
}

getT = function(n,k,x=1:n) {
  T = x
  if (k==0) T = T[-1]
  else if (k%%2==1) T = T[-c(1:((k+1)/2),(n-(k-1)/2):n)]
  else if (k%%2==0) T = T[-c(1:((k+2)/2),(n-(k-2)/2):n)]
  return(T)
}

getG = function(n,k,x=1:n) {
  T = getT(n,k,x)
  G = matrix(0,n,n)
  G[,1] = rep(1,n)
  for (j in 2:n) {
    if (j<=k+1) {
      G[,j] = x^(j-1)
    }
    else {
      G[,j] = pmax(x-T[j-k-1],0)^k
    }
  }
  return(G)
}


# Correction because getH() for trendfiltering above k>0 is actually wrong;
# This is from Ryan's trendfilter paper Lemma 2
getH.trendfilter = function(n,k,x=1:n){

  # get k'th order cumulative sum
  s = rep(1,n)
  if(k>=1){
    for(ik in 1:k)   s = cumsum(s)
  }

  # Fill in H matrix
  H = matrix(NA,nrow=n,ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      if(j<=(k+1)){
        H[i,j] =  (i/n)^(j-1)
      } else {
        if(i <= j-1) {
          H[i,j] = 0
        } else {
          H[i,j] = s[i-j+1]*factorial(k)/(n^k)
        }
      }
    }
  }

  return(H)
}


getHslow = function(n,k) {
  D = matrix(0,n,n)
  D[1,1] = 1
  for (i in 2:(k+1)) D[i,] = getDtf(n,i-2)[1,]
  D[(k+2):n,] = getDtf(n,k)
  return(round(solve(D),3))
}

getHscaled = function(n,k,x=1:n) {
  H = getH(n,k,x)
  return(H*factorial(k))
}

compare = function(a,aa,lam) {
  f = predict(a,lambda=lam)$fit
  if ("trendfilter" %in% class(aa)) ff = coef(aa,lambda=lam)$beta
  else ff = predict(aa,lambda=lam)$fit
  n = length(a$y)
  plot(1:n,a$y)
  lines(1:n,f,col="blue")
  lines(1:n,ff,col="red")
  legend("bottomleft",col=c("blue","red"),
         legend=c("LARS","TF"),lty=1)
  return(max(abs(f-ff)))
}

kspline = function(k=0, n=100, seed=0, numknots=5, knots=NULL, weights=NULL) {
  if (is.null(knots)) knots=sample(1:n,numknots)
  else numknots = length(knots)

  if (is.null(weights))
    weights=sample(c(-1,1),numknots+k+1,replace=TRUE)*rnorm(numknots+k+1,2,1)

  beta = rep(0,n)
  for (i in 0:k) {
    beta = beta + weights[i+1]*(1:n)^k/n^k
  }
  for (j in 1:numknots) {
    x = 1:n-knots[j]
    x = x^k*(x>0)
    beta = beta + weights[j+k+1]*x/max(x)
  }

  beta = 10*(beta-min(beta))/(max(beta)-min(beta))
  y = beta + rnorm(n,0,1)

  return(list(beta=beta,y=y))
}

pieconst = function(n=100, seed=0, numknots=5) {
  set.seed(seed)

  d = floor(n/numknots)
  beta = rep(sample(1:10,numknots),
    times=c(rep(d,numknots-1),n-(numknots-1)*d))
  beta = 10*(beta-min(beta))/(max(beta)-min(beta))

  y = beta + rnorm(n,sd=1)
  return(list(beta=beta,y=y))
}

pielin = function(n=100, seed=0) {
  set.seed(seed)
  n = 100
  knots = matrix(0,6,2)
  knots[,1] = c(1,20,45,60,85,100)
  knots[,2] = c(1,2,6,8,5,6)
  beta = rep(0,100)
  for (i in 1:(nrow(knots)-1)) {
    for (j in knots[i,1]:(knots[i+1,1])) {
      beta[j] = (knots[i+1,1]-j)/(knots[i+1,1]-knots[i,1])*knots[i,2] +
        (j-knots[i,1])/(knots[i+1,1]-knots[i,1])*knots[i+1,2]
    }
  }

  beta = 10*(beta-min(beta))/(max(beta)-min(beta))
  y = beta+rnorm(n,0,1)

  return(list(beta=beta,y=y))
}

pielinx = function(x) {
  knots = matrix(0,6,2)
  knots[,1] = c(0,20,45,60,85,100)/100
  knots[,2] = c(1,2,6,8,5,6)
  ii = max(which(x>=knots[,1]))
  if (ii==nrow(knots)) {
    return(knots[nrow(knots),2])
  }
  else {
    a = (x-knots[ii,1])/(knots[ii+1,1]-knots[ii,1])
    return(knots[ii,2]*(1-a) + knots[ii+1,2]*a)
  }
}


piequad = function(n=100, seed=0) {
  set.seed(seed)
  knots = matrix(0,4,2)
  knots[,1] = c(1,33,60,100)
  knots[,2] = c(8,6,2,4)
  beta = rep(0,n)
  endval = 0
  for (i in 1:(nrow(knots)-1)) {
    sgn = (-1)^(i+1)
    mid = (knots[i,1]+knots[i+1,1])/2
    dif = knots[i+1,1]-knots[i,1]

    j = knots[i,1]
    intcp = endval - (sgn*n/5*(j-mid)^2/dif^2 + (knots[i+1,1]-j)/dif*knots[i,2] +
      (j-knots[i,1])/dif*knots[i+1,2])

    for (j in knots[i,1]:(knots[i+1,1])) {
      beta[j] = intcp + sgn*n/5*(j-mid)^2/dif^2 + (knots[i+1,1]-j)/dif*knots[i,2] +
        (j-knots[i,1])/dif*knots[i+1,2]
    }
    endval = beta[j]
  }

  beta = 10*(beta-min(beta))/(max(beta)-min(beta))
  y = beta + rnorm(n,0,1)

  return(list(beta=beta,y=y))
}

piecub = function(n=100, seed=0) {
  set.seed(seed)
  n=100
  beta=rep(0,100)
  beta[1:40]=(1:40-20)^3
  beta[40:50]=-60*(40:50-50)^2 + 60*100+20^3
  beta[50:70]=-20*(50:70-50)^2 + 60*100+20^3
  beta[70:100]=-1/6*(70:100-110)^3 + -1/6*40^3 + 6000
  beta=-beta

  beta=10*(beta-min(beta))/(max(beta)-min(beta))
  y=beta+rnorm(n,0,1)

  return(list(beta=beta,y=y))
}

piecubx = function(x) {
  x = x*100
  if (x<=40) return((x-20)^3)
  if (x<=50) return(-60*(x-50)^2 + 60*100+20^3)
  if (x<=70) return(-30*(x-50)^2 + 60*100+20^3)
  if (x<=100) return(-1/6*(x-110)^3 + -1/6*40^3 + 6000)
}

smoothwiggly.fun = function(x) {
  f = function(a,b,c,x) return(a*(x-b)^2+c)
  fp = function(a,b,c,x) return(2*a*(x-b))

  a=-1; b=1/4; c=1;
  if (x<=1/3) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=1/3;
  a=1; b=xx-fp(aa,bb,cc,xx)/(2*a); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=2/3) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=2/3;
  b=0.7; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.775) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.775;
  b=0.8; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.825) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.825;
  b=0.85; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.875) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.875;
  b=0.9; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.925) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.925;
  b=0.95; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.975) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.975;
  b=1; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  return(f(a,b,c,x))
}

smoothwiggly = function(n=100, x=1:n/n) {
  u = rep(0,n)
  for (i in 1:n) u[i]=smoothwiggly.fun(x[i])
  return(u)
}

maxlam2 = function(y,k) {
  n = length(y)
  D = getDtfSparse(n,k)
  x = qr(t(D))
  u = backsolveSparse(x,y)
  return(list(u=u,maxlam=max(abs(u))))
}

## maxlam3 = function(y,k) {
##   n = length(y)
##   D = getDtfSparse(n,k)
##   x = qr(crossprod(t(D)))
##   u = backsolveSparse(x,D%*%y)
##   return(list(u=u,maxlam=max(abs(u))))
## }

backsolveSparse = function(QR, b) {
  #R = qrR(QR)
  #x = solve(R, qr.qty(QR,b)[Seq(1,nrow(R))])
  #return(as.numeric(x))
  R = qr.R(QR)
  x = solve(R, qr.qty(QR,b)[Seq(1,nrow(R))])
  if (length(QR@q)==0) return(x)
  else return(x[Order(QR@q+1)])
}

Order = function(x) {
  n = length(x)
  o = numeric(n)
  o[x] = Seq(1,n)
  return(o)
}

crit = function(y,k,lambda,beta) {
  return(0.5*sum((y-beta)^2) + lambda*sum(abs(diff(beta,differences=k+1))))
}

bx = function(x,lam) {
  y = x
  y[x>lam] = lam
  y[x< -lam] = -lam
  return(y)
}

st = function(x,s,r) {
  z = rep(0,length(x))
  z[x>s] = x[x>s]-s
  z[x< -s] = x[x< -s]+s
  # leave the first r coordinates untouched
  z[1:r] = x[1:r]
  return(z)
}

# Make sure that no two columms of A are the same
# (this works with probability one).

checkcols <- function(A) {
  b = rnorm(nrow(A))
  a = sort(t(A)%*%b)
  if (any(diff(a)==0)) return(TRUE)
  return(FALSE)
}

##' Helper functions for fl1d (by Justin)
LS = function(y, X, out = c("fitted_y", "fitted_coeff")){
  ## n-length vector y
  ## n by p matrix X  (n >= p)
  if( out == "fitted_y" ){
    Projection <- X %*% solve( t(X) %*% X ) %*% t(X)
    yhat <- Projection %*% y
    return(yhat)
  }
  if( out == "fitted_coeff"){
    Beta_hat <- solve( t(X) %*% X ) %*% t(X)
    return(Beta_hat)
  }
}


dual1d_Dmat = function(m){
  D = matrix(0, nrow = m-1, ncol = m)
  for(ii in 1:(m-1)){
    D[ii,ii] = -1
    D[ii,ii+1] = 1
  }
  return(D)
}

# Makes a D matrix for a matrix that is stacked as rows
graph2D_Dmat = function(m){

    m0 = sqrt(m)
    stopifnot( m0 == round(m0) )
    D = matrix(0, nrow = 2*(m0-1)*m0, ncol = m)

    # connect all within chunks
    for(jj in 1:m0){
      chunk = 1:m0 + (jj-1)*m0
      D[(m0-1)*(jj-1) + 1:(m0-1), chunk] = dual1d_Dmat(m0)
    }
    sofar = m0*(m0-1)
    # connect between all adjacent chunks
    for(jj in 1:(m-m0)){
      newrow = rep(0,m)
      newrow[jj]    = -1
      newrow[jj+m0] =  1
      D[sofar+jj,] = newrow
    }
  return(D)
}

##' general function that makes various D matrices
##' @export
makeDmat = function(m, type = c("tf","2d","graph"), order=0){

    type = match.arg(type)
    D = NA
    if(type == "tf"){
        ## Multiplies D matrices to form k'th order tf matrix
        D = dual1d_Dmat(m)
        if(order>=1){
            for(jj in 1:order){
                D = dual1d_Dmat(m-jj) %*% D
            }
        }
    } else if ("2d") {
        D = graph2D_Dmat(m)
    } else {
        stop("Not coded yet!")
    }
    return(D)
}

# helper function, to make a vector input into a row matrix
asrowmat = function(obj){
  objmat = as.matrix(obj)
  if(ncol(objmat)==1 & ncol(objmat) < nrow(objmat)){
    objmat = t(objmat)
  }
  return(objmat)
}

# and likewise, column matrix
ascolmat = function(obj){
  objmat = as.matrix(obj)
  if(nrow(objmat)==1 & nrow(objmat) < ncol(objmat)){
    objmat = t(objmat)
  }
  return(objmat)
}




mod = function(num,mod){
  c = num %% mod
  if(c == 0){
    return  (mod)
  }
  else{
    return (c)
  }
}

## from github repository statsmath/genlasso/R/svdsolve.R, changing default rtol
svdsolve <- function(A,b, rtol=1e-7) {
  s = svd(A)
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v%*%(di*(t(s$u)%*%b)),q=sum(ii),s=s))
}


#' creates vector of length n whose j'th element is 1, otherwise 0. Not currently in use.
evec = function(n,j){
  v <- rep(0,n)
  v[j] <- 1
  return(v)
}

# Transform the fused lasso problem into a lasso problem, using (11) of Genlasso path paper
getlars = function(y,k,numsteps){
  n = length(y)
  H = getH(n,k)
  H1 = H[,1:(k+1)]
  H2 = H[,(k+2):n]
  Py = H1%*%solve(crossprod(H1),t(H1)%*%y0)
  yy = y0 - Py
  PH = H1%*% solve(crossprod(H1),t(H1)%*%H2)
  HH = H2 - PH

  # extract gamma from lasso path
  out = lar(yy,HH,maxsteps=numsteps,verbose=FALSE)
  return(out)
}




#' produces primal solution from dual solution on a point in path
getprimal_from_dual = function(u, y, type = "fl1d"){
  if(type == "fl1d"){
      D = dual1d_Dmat(length(y))
  } else {
      stop("type not coded yet")
  }
  return( y - t(D)%*%u )
}

#' produces primal matrix from dual solution matrix
#' Example usage: getprimal(umat = fl1d(..)$u)
getprimalmat = function(umat, y, type="fl1d"){
  if(type == "fl1d"){
    D = dual1d_Dmat(length(y))
  } else {
    stop("type not coded yet")
  }
  bmat = apply(X = umat, MARGIN = 2, function(u){ return(getprimal(u,y, type)) })# 2 to apply over columns
  return(bmat)
}

# Make residuals from the differences in the null space projection
makeresid = function(D,prevhits,thishit,n){

      proj = function(mymat){ return(mymat %*% solve(t(mymat)%*%mymat, t(mymat)))}

      # Record t(D_{-B}) for previous and current models
      tD.curr = t(D[,-c(prevhits,thishit)])
      tD.prev = t(D[,-prevhits])
      tDb.curr = svd(tD.curr)$u[,1:rankMatrix(tD.curr)]
      tDb.prev = svd(tD.prev)$u[,1:rankMatrix(tD.prev)]

      # form orth projections to row spaces
      p.curr = diag(1,n) - proj(tDb.curr)
      p.prev = diag(1,n) - proj(tDb.prev)
      p.resid = (pprev - pcurr)

      # extract the first projection basis
      myresid = svd(p.resid)$u[,1]
      return(myresid)
}

# Make residuals from projection, from equivalent regression problem
# prevhits and thishit should have no overlap
makeresid.tf = function(prevhits, thishit,n, k=0){

    # Get old model (H1) and new variable being added (H2)
    H = getH(n,k)
    nullmodel.basis.inds = 1:(k+1)
    H1 = (if(is.na(prevhits[1])) H[,nullmodel.basis.inds] else H[,c(nullmodel.basis.inds,prevhits+1)])
    H2 = H[,(thishit+ 1)]

    # Obtain residual from regressing additional variable to small(old) model
    resid = resid(lm(H2 ~ H1-1))

  return(resid)
}


# Makes the residual vector from projecting a new group
# membership onto the old ones
# |new.membership|: index set of members of a newly added group.
# |old.memberships|: matrix whose columns are partial-1-vectors containing
#                    the existing graph cluster memberships.
# For example:
#new.membership = c(1,1,0,0,0,0)
#old.memberships = cbind(c(1,1,1,1,1,1),c(1,1,1,1,0,0))
makeresid.graph = function(new.membership, old.memberships, n){

  # make sure to check if new.index is a subset of old.index
  stopifnot(all(as.numeric(old.memberships - new.membership) %in% c(0,1)))

  # make residual vector by same argument of projecting new membership on old
  residual = resid(lm(new.membership ~ old.memberships))
  return(residual)
}


makeresid = makeresid.tf




##' This function creates a color scale for use with the image() function. Input
##' parameters should be consistent with those used in the corresponding image
##' plot. The "axis.pos" argument defines the side of the axis. The "add.axis"
##' argument defines whether the axis is added (default: TRUE)or not (FALSE).
##' Taken from http://www.r-bloggers.com/new-version-of-image-scale-function/
image.scale <- function(z, zlim, col = heat.colors(12),
breaks, axis.pos=1, add.axis=TRUE, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
 if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
 plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
 for(i in seq(poly)){
  if(axis.pos %in% c(1,3)){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(axis.pos %in% c(2,4)){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
 box()
 if(add.axis) {axis(axis.pos)}
}



##' Function to show top n sized objects and their sizes
show.glutton <- function(envir=.GlobalEnv, how.many = 5){
  objs=ls(envir=envir)
  obj.sizes.byte = sapply(objs,function(obj){object.size(get(obj,envir=envir))})
  obj.sizes.Mb = sapply(objs,function(obj){format(object.size(get(obj,envir=envir)),"Mb")})
  obj.sizes = head(obj.sizes.Mb[  order(obj.sizes.byte,decreasing=T)],how.many)
  print(obj.sizes)
}



##' Function to show objects larger than 10 Mb in memory
see.big.objs <- function(szcut=10, which.envir = c("global", "local"), verbose=F){
  which.envir = match.arg(which.envir)
  if(which.envir=="global"){
    envir = .GlobalEnv
   } else if (which.envir == "local"){
    envir = sys.frame()
   } else {
    stop("Enviroment not properly specified!")
   }

  large.sz = sapply(ls(envir=envir), function(objname) { sz = gsub('.{3}$', '', format(object.size(get(objname,envir=envir)),"Mb"));  if(as.numeric(sz)>szcut){return(sz)} else {return(NA)}  })
  large.sz = large.sz[!is.na(large.sz)]
  if(verbose) print(large.sz)
  return(large.sz)
}



##' Function to take mypath$action from a path object and returns a list of the
##' states after each step in algorithm, starting with NA as the first state.
##' @export
get.states = function(action.obj, path.obj=NULL){

    ## Extract action from pathobj
    if(!is.null(path.obj)){
        action.obj = path.boj$action
    }

    ## Obtain a list of actions at each step.
    actionlists = sapply(1:length(action.obj),function(ii) action.obj[1:ii] )

    ## Helper function to extract the final state after running through (a list of) actions
    get.final.state = function(actionlist){
        if(length(actionlist)==1) return(actionlist)
        to.delete = c()
        for(abs.coord in unique(abs(actionlist))){
            all.coord.actions = which(abs(actionlist)==abs.coord)
            if( length(all.coord.actions) > 1 ){
                if( length(all.coord.actions) %%2 ==1){
                    to.delete = c(to.delete, all.coord.actions[1:(length(all.coord.actions)-1)])
                } else {
                    to.delete = c(to.delete, all.coord.actions)
                }
            }
        }
        if(!is.null(to.delete)){
            return(actionlist[-to.delete])
        } else {
            return(actionlist)
        }
    }

    ## Extract them
    states = lapply(actionlists, get.final.state)
    states = c(NA,states)
    return(states)
}

##' Helper function to obtain the effective design X matrix, by breaking X at
##' the fused lasso regression breakpoints.
get.augmented.X = function(X,breaks,return.groups = F,TT,J,group.inds){
    find.which.group = function(hit) which(sapply(1:J, function(jj) hit %in% group.inds[[jj]]))
    which.group.breaks = sapply(breaks, find.which.group)
    breaks.by.group = list()
    for(groupnum in 1:J){
      thisgroup.breaks = breaks[which.group.breaks==groupnum]
      thisgroup.breaks = thisgroup.breaks - TT*(groupnum-1) + (groupnum-1)
      breaks.by.group[[groupnum]] = thisgroup.breaks
    }

    # Function to break var into #|breaks| variables broken at |breaks|
    brk.var = function(var,breaks){
      augmented.breaks =c(0,sort(breaks),length(var))
      all.subvecs =  sapply(1:(length(augmented.breaks)-1),
                     function(ii){inds = (augmented.breaks[ii]+1):(augmented.breaks[ii+1])
                                  newvar = rep(0,length(var))
                                  newvar[inds] = var[inds]
                                  return(newvar)
                                 })
      return(all.subvecs)
    }

    X.augmented = do.call(cbind,
                          lapply(1:length(breaks.by.group),
                                  function(ii) brk.var(X[,ii], breaks.by.group[[ii]])))

     if(return.groups){ return(breaks.by.group) } else { return(X.augmented) }

  }

##' Returns sequence of BIC of fused lasso regression model (general regressor
##' matrix X) The variable names 00.orig are to emphasize that while the f0 may
##' be from solving a ridge-penalty-added problem, what should be provided are
##' the original, non-augmented versions
getbic.regression = function(y0.orig, f0, sigma, maxsteps,X.orig, ginvX.orig, D.orig,rtol=1E-7){

  ## Setup empty IC vector / problem dimensions / pseudoinverse.
  bic = numeric(maxsteps+1)
  TT=length(y0.orig) # largest'time'
  J = ncol(X.orig)   # number of variables in the regression
  n = length(y0.orig)
  if(is.null(ginvX.orig)) ginvX.orig = ginv(X.orig)

  # Obtain each path's state at each step
  states = get.states(f0$action)

  # Collect BIC at each step 0 ~ (maxsteps-1)
  for(ii in 1:(maxsteps+1)){

######### OLD WAY:  Find the basis of the row space of D_{-B}  ##################################################################3
#    current.state = states[[ii]]
#    D.interior = (if(ii==1) D else D[-current.state,])
#    proj.orth.row.space = diag(rep(1,ncol(D.interior))) - t(D.interior)%*%solve(D.interior %*% t(D.interior), D.interior)
#    tH = t(D.interior %*% ginv(X.orig)) # basis for null space of D_(-B)
#    y0.fitted = as.numeric(resid(lm(y0~tH-1)))
#    y0.fitted2 = proj.orth.row.space %*% y0

######## NEW WAY:  Find the basis of the row space of D_{-B}  ##################################################################3
    curr.state = states[[ii]]

    # Function to obtain the "augmented" X matrix, by breaking at the fused lasso breakpoints
    get.augmented.X = function(X.orig,breaks,return.groups = F){
      group.inds = lapply(1:J, function(ii){(TT*(ii-1)+1):(TT*(ii)-1) - (ii-1) })
      find.which.group = function(hit) which(sapply(1:J, function(jj) hit %in% group.inds[[jj]]))
      which.group.breaks = sapply(breaks, find.which.group)
      breaks.by.group = list()
      for(groupnum in 1:J){
        thisgroup.breaks = breaks[which.group.breaks==groupnum]
        thisgroup.breaks = thisgroup.breaks - TT*(groupnum-1) + (groupnum-1)
        breaks.by.group[[groupnum]] = thisgroup.breaks
      }

      # Function to break var into #|breaks| variables broken at |breaks|
      brk.var = function(var,breaks){
        augmented.breaks =c(0,sort(breaks),length(var))
        all.subvecs =  sapply(1:(length(augmented.breaks)-1),
                       function(ii){inds = (augmented.breaks[ii]+1):(augmented.breaks[ii+1])
                                    newvar = rep(0,length(var))
                                    newvar[inds] = var[inds]
                                    return(newvar)
                                   })
        return(all.subvecs)
      }

      X.augmented = do.call(cbind,
                            lapply(1:length(breaks.by.group),
                                   function(ii) brk.var(X.orig[,ii], breaks.by.group[[ii]])))
      if(return.groups){ return(breaks.by.group) } else { return(X.augmented) }
    }

    # Augment the design matrix at their breaks
    X.aug.curr  = get.augmented.X(X.orig, curr.state )
    sv = svd(X.aug.curr)
    Xrank = invisible(rankMatrix(X.aug.curr))
    b = sv$u[,1:Xrank]
    # Get the basis vector of the residual linear subspace
    projection = function(mat){
      mat %*% solve(t(mat) %*% mat, t(mat))
    }
    Proj.curr = projection(b) #X.aug.curr %*% solve(t(X.aug.curr)%*% X.aug.curr, t(X.aug.curr))
    #y0.fitted  = Proj.curr %*% (y0.orig)
    y0.fitted = as.numeric(fitted(lm(y0.orig ~ X.aug.curr-1)))
    beta.fitted = as.numeric(coef(lm(y0.orig ~ X.aug.curr-1)))
    mygroups = get.augmented.X(X.orig,curr.state,T)


#    # sanity check
#    beta.fitted.long = rep(NA,TT*J)
#    beta.fitted.long[(1):(TT)] = beta.fitted[1]
#    beta.fitted.long[(TT+1):(TT+mygroups[[2]][1])] = beta.fitted[2]
#    beta.fitted.long[(TT+mygroups[[2]][1]+1):(2*TT)] = beta.fitted[3]
#    y0.fitted2 = X.augmented %*% beta.fitted.long
#    plot(y0.fitted,col='red');lines(y0.fitted2,col='red');points(y0,pch=16)

    ## Calculate bic
      RSS = sum((y0.orig - y0.fitted)^2)
      ## cat("RSS is", RSS, fill=T)
      getdf = function(D,rtol=1e-7){
          mysv = svd(diag(rep(1,ncol(D))) - projection(t(D)))
          return(sum(mysv$d>rtol))
      }
      rowind = (if(all(is.na(curr.state))) (1:nrow(D.orig)) else (1:nrow(D.orig))[-curr.state[!is.na(curr.state)]])
      df = getdf(D.orig[rowind,])#sum(sv$d > rtol)#invisible(rankMatrix(Proj.curr))
      ## cat("df is", df, fill=T)
      complexity.multiplier = sigma^2 * log(TT*J+TT)
      complexity.penalty = df * complexity.multiplier
      ## cat("multipl is ", complexity.multiplier, fill=T)
    bic[ii] <- RSS + complexity.penalty
  }
  return(bic)
}


##' Calculates various model-related quantities for _any_ D for the signal
##' approximator case (IC, penalty, RSS, residual) \code{df.fun()} is a function
##' that returns the degrees of freedom of fit of %yhat = Proj_{null(D_B)} y%.
##' Since there is no reason to /want/ to use a different \code{D} matrix than
##' what was used in the path object \code{f0}, we don't consider
##' @param f0 Object of class |path|.
##' @param consec How many consecutive rises in IC? Defaults to 2
##' @param sigma Standard deviation that the i.i.d. (Gaussian) data was
##'     generated from.
##' @param maxsteps Maximum number of path steps to consider, in calculating
##'     degrees of freedom. Defaults to the total number of steps made in the
##'     path (in \code{f0}), in the first place!
##' @param stoprule Out of BIC, EBIC, and AIC, which criterion you want to
##'     consider.
##' @param ebic.fac Correction factor in EBIC by Chen & Chen (2008)
##' @param verbose Whether to be loud while doing stuff.
##' @export
get.modelinfo = function(obj, consec=2, sigma, maxsteps=length(obj$action),
                         stoprule = c('bic','ebic', 'aic'), ebic.fac=.5,
                         verbose=F){

    ## Basic checks
    stoprule = match.arg(stoprule)

    ## Load things
    D = obj$D
    y0 = obj$y
    n = length(y0)

    ## Create empty objects
    ic = RSS = pen = df = rep(NA, maxsteps)
    resids = matrix(NA,nrow=n,ncol=maxsteps)
    actiondirs = c(NA, sign(obj$action)[1:(maxsteps)])  # s* : change in model size

    ## Obtain sequence of path's state at the end of each step
    states = get.states(obj$action)

    ## Collect BIC at each step 0 ~ (maxsteps-1)
    for(ii in 1:maxsteps){
        if(verbose){
            cat('step', ii, '\n')
        }
        ## Method 1: Form proj null(D_{-B}) by orth proj onto row(D_{-B}) = col(t(D_{-B})) ~= tD
        thishits = states[[ii]]
        tD = cbind(if(all(is.na(thishits))) t(D) else t(D)[,-(thishits)])

        proj = function(mymat){ return(mymat %*% solve(t(mymat)%*%mymat, t(mymat)))}
        rr = rankMatrix(tD)
        tDb = svd(tD)$u[,1:rr]
        curr.proj = proj(tDb)
        y0.fitted = (diag(1,n) - curr.proj) %*% y0

        ##      if(ii>1) {
        ##        prevhits = states[[ii-1]]
        ##        tD.prev = cbind(if(all(is.na(prevhits))) t(D) else t(D)[,-(prevhits)])
        ##        rr.prev = rankMatrix(tD.prev)
        ##        tDb.prev = svd(tD.prev)$u[,1:rr.prev]
        ##        prev.proj = proj(tDb.prev)
        ##        y0.fitted.prev = (diag(1,n) - prev.proj) %*% y0
        ##      }
        ##
        ##      par(mfrow=c(2,1))
        ##      plot(y0.fitted,ylim = range(y0),type='o'); lines(y0.fitted.prev,col='red');
        ##      plot(logb(bic,base=100),ylim=c(-1,1)); abline(v=ii,col='red')
                                        #
        ## Method 2: Form null projection onto D_{-B} directly, from svd
        ##      thishits = states[[ii]]
        ##      D.curr = if(all(is.na(thishits))) D else D[-(thishits+1),]
        ##      rD.curr = rankMatrix(D.curr)
        ##      if(rD.curr == ncol(D.curr)){
        ##        null.proj = Matrix(0,nrow = length(y0),ncol = length(y0))
        ##      }
        ##      else {
        ##        null.curr = svd( D.curr, nv = ncol(D.curr))$v[,(rD.curr+1):ncol(D.curr)]
        ##        null.proj = null.curr %*% t(null.curr)
        ##      }
        ##      y0.fitted2 = null.proj %*% y0

        myRSS = sum( (y0 - y0.fitted)^2 )
        mydf  = n-rr
        if(ii==1) prev.df = mydf
        mypen = (if(stoprule=='bic'){
                     (sigma^2) * mydf * log(n)
                 } else if (stoprule=='ebic'){
                     (sigma^2) * mydf * log(n) + 2*(1-ebic.fac) *
                         log(choose(n,mydf)) * sigma^2
                 } else if (stoprule=='aic'){
                     (sigma^2) * mydf * 2
                 } else {
                     stop(paste(type, "not coded yet!"))
                 })

                                        # If there is an addition to the boundary set and an addition to the
        if(ii>1  & mydf != prev.df){#& actiondirs[ii] == +1
                                        #        rankMatrix(curr.proj - prev.proj)
            ptol = 1E-10
            if(all(abs(as.numeric(curr.proj - prev.proj)) < ptol)){
                myresid = rep(NA,n)
            } else {
                myresid = svd(curr.proj - prev.proj)$u[,1]
                myresid = myresid / sqrt(sum((myresid)^2))
            }
        } else {myresid = rep(NA,n)}

        ## Store RSS + p(df)
        RSS[ii] <- myRSS
        pen[ii] <- mypen
        ic[ii] <- myRSS + mypen
        resids[,ii] <- myresid
        df[ii] <- mydf

        ## Update previous information by current information.
        prev.proj = curr.proj
        prev.df = mydf
    }

    ## Record things at primal-changing knots (useful for the graph case)
    ic.primal = rep(NA,maxsteps)
    ic.primal[1] = ic[1]
    df.before = 1
    for(ii in 2:length(df)){
        if(!(df.before == df[ii])){
            ic.primal[ii] = ic[ii]
        }
        df.before = df[ii]
    }
    knots.primal  = which(!is.na(ic.primal))
    resids.primal = resids[knots.primal]
    ic.primal     = ic.primal[knots.primal]
    actiondirs.primal = actiondirs[knots.primal]

    ## Issue warning if BIC hasn't stopped, or stoptime is zero.
    stop.time = which.rise(ic.primal,consec=consec) - 1
    stop.time = pmin(stop.time,n-consec-1)
    stopped=TRUE
    if(!(stop.time+consec < maxsteps)){
        warning(paste('IC rule using', consec, 'rises hasnt stopped!'))
        stopped = FALSE
    }
    stop.time = pmin(stop.time, n-consec-1) ## can't remember why I did this.
    if(stop.time==0) warning('Stoptime is zero!')

    return(list(RSS=RSS, pen=pen, ic=ic, resids=resids,
                knots.primal = knots.primal, resids.primal = resids.primal,
                ic.primal = ic.primal, actiondirs.primal = actiondirs.primal,
                stopped = stopped))
}



## Calculates BIC for _any_ D for the signal approximator case df.fun is a
## function that returns the degrees of freedom of fit of yhat =
## Proj_{null(D_B)} y
get.ic = function(obj, y0, sigma, maxsteps, D, type = c('bic','ebic', 'aic'), ebic.fac=.5){

  # Setup empty IC vector and regressor matrix H, for k'th order trend filtering problems
  ic = RSS = pen = numeric(maxsteps+1)
  n = length(y0)

  # Obtain each path's state at each step
  states = get.states(obj$action)

  # Collect BIC at each step 0 ~ (maxsteps-1)
  for(ii in 1:(maxsteps+1)){

    # Method 1: Form proj null(D_{-B}) by orth proj onto row(D_{-B}) = col(t(D_{-B})) ~= tD
      thishits = states[[ii]]
      tD = cbind(if(all(is.na(thishits))) t(D) else t(D)[,-(thishits)])
      proj = function(mymat){ return(mymat %*% solve(t(mymat)%*%mymat, t(mymat)))}
      rr = rankMatrix(tD)
      tDb = svd(tD)$u[,1:rr]
      y0.fitted = (diag(1,n) - proj(tDb)) %*% y0

    # Method 2: Form null projection onto D_{-B} directly, from svd
#      thishits = states[[ii]]
#      D.curr = if(all(is.na(thishits))) D else D[-(thishits+1),]
#      rD.curr = rankMatrix(D.curr)
#      if(rD.curr == ncol(D.curr)){
#        null.proj = Matrix(0,nrow = length(y0),ncol = length(y0))
#      }
#      else {
#        null.curr = svd( D.curr, nv = ncol(D.curr))$v[,(rD.curr+1):ncol(D.curr)]
#        null.proj = null.curr %*% t(null.curr)
#      }
#      y0.fitted2 = null.proj %*% y0

    myRSS = sum( (y0 - y0.fitted)^2 )

    mydf = n-rr

    if(type=='bic')  mypen = (sigma^2) * mydf * log(n)
    if(type=='ebic') mypen = (sigma^2) * mydf * log(n) +
                              2*(1-ebic.fac) * log(choose(n,mydf)) * sigma^2
    if(type=='aic')  mypen = (sigma^2) * mydf * 2

    myresid = #residual from adding the latest variable to the existing variable

    # Store RSS + p(df)
    RSS[ii] <- myRSS
    pen[ii] <- mypen
    ic[ii] <- myRSS + mypen
  }
  return(list(RSS=RSS,pen=pen,ic=ic))()
}

# Returns a sequence of +1 and -1 for sequential incr and decrements in a vector
# assume first step always dips; _almost_ always true
getorder = function(bic){
  return(c(NA,sign(bic[2:(length(bic))] - bic[1:(length(bic)-1)])))
}


##' Locates first point of .(consec) rises in IC path
##' @export
which.rise = function(icvec, consec = 2, direction=c("forward","backward")){

  direction = match.arg(direction)
  if(direction != "forward") stop("That direction IC selection is not coded yet.")
  if(length(icvec) < consec+1) stop("Not enough steps to do forward sequential BIC/AIC!")

  ind = 1
  done = FALSE
  while(ind < (length(icvec)-consec+1) ){
    ictol = 1E-10
    if( all(icvec[(ind+1):(ind+consec)] > icvec[(ind):(ind+consec-1)] + ictol )) break
    ind = ind+1
  }
  return(pmin(pmax(ind,1),length(icvec)))
}

##' Takes in adjacency matrix (either produced from igraph, or a matrix with
##'                            zero entries in non-adjacent edges and 1 in
##'                            adjacent edges) and produces a (sparse) D matrix
##'                            for usage in the generalized lasso
##' @import Matrix
getDmat.from.adjmat = function(adjmat, sparseMatrix=FALSE){
  n = ncol(adjmat)
  if(sparseMatrix){
    Dmat = Matrix(0,nrow=n^2,ncol=n, sparse=TRUE)
  } else {
    Dmat = matrix(0,nrow=n^2,ncol=n)
  }
  count = 1
  for(jj in 1:n){
    for(kk in jj:n){
        if(adjmat[jj,kk]==1){ # there was an error !=1
          Dmat[count,jj] = 1
          Dmat[count,kk] = -1
          count = count+1
      }
    }
  }
  Dmat = Dmat[1:(count-1),]
  return(Dmat)
}


##' Takes in a D matrix used in the generalized lasso and creates an adjacency
##' matrix (a matrix with zero entries in non-adjacent edges and 1 in adjacent
##' edges).
getadjmat.from.Dmat = function(Dmat){

  Dmat = rbind(Dmat)
  n = ncol(Dmat)
  adjmat = matrix(0,n,n)

  for(jj in 1:nrow(Dmat)){
    inds = which(Dmat[jj,]!=0)
    adjmat[inds[1],inds[2]] = 1
    adjmat[inds[2],inds[1]] = 1
  }
  return(adjmat)
}


## takes in actions from a dualpathsvd(2) object and returns the boundary set at that step
getB.from.actions = function(actions){
  B = actions
  to.delete = c()
  for(ww in 1:length(actions)){
    if(actions[ww] < 0){
      to.delete = c(to.delete, which(actions == abs(actions[ww])))
      to.delete = c(to.delete, ww)
    }
  }
  if(length(to.delete)!=0){
    B = B[-to.delete]
  }
  return(B)
}

##'  Wrapper for getmodelinfo.graph() to get bic scores only Returns BIC for
##'  steps 0~maxsteps
getbic.graph = function(obj,y0, sigma, maxsteps, Dmat, stoptime=F, dual = T){
  a = getmodelinfo.graph(obj,y0, sigma, maxsteps, Dmat, stoptime=F)
  if(dual){ return(a$bic) } else { return(a$bic.primal)}
}

##'  Function to create things related to model selection, for the graph case
##'  Input : path object, y0, initial graph, maximum steps to take, Dmat,
##'  cluster size, cluster # Output: 6 p values for each possible segment test.
getmodelinfo.graph = function(obj,y0, sigma, maxsteps, Dmat, stoptime=F){

  mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat), mode="undirected")
  n = length(y0)
  oldmodel = cbind(rep(1,n))
  models = list()
  df = c()
  residuals = list()

  # Get initial group membership + store null model
  prev.groups = list()
  for(member in unique(igraph::clusters(mygraph)$membership))  prev.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
  prev.nclust = length(prev.groups)
  models[[1]] = oldmodel
  df[1] = 1
  residuals[[1]] = NA

  for(step in 1:maxsteps){
      ## Get current Dmat
      Dmat.curr = Dmat[-getB.from.actions(obj$action[1:step]),]

      ## Get current Graph
      mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat.curr), mode="undirected")
      if(nrow(rbind(Dmat.curr))==0) break
      num.connected.components = length(unique(igraph::clusters(mygraph)$membership))

      ## Store degrees of freedom of this step
      df[step+1] = num.connected.components

      ## Keep track of connected components
      curr.groups = list()
      for(member in unique(igraph::clusters(mygraph)$membership)){
          curr.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
      }
      curr.nclust = length(curr.groups)

      newbasis = rep(0,length(y0))

      ## If cluster was created
      if(curr.nclust == prev.nclust + 1){
          if(verbose==TRUE) { print("cluster created!") }
          matched.curr.group.ind = c()
                                        # Which cluster was created? Scan all groups from prev and curr; what has changed?
          for(jj in 1:prev.nclust){
              for(kk in 1:curr.nclust){
                  if( setequal(prev.groups[[jj]], curr.groups[[kk]])){
                      ## if _any_ of the previous groups match with current groups exactly, assign zero
                      newbasis[prev.groups[[jj]]] = 0
                      matched.curr.group.ind = c(matched.curr.group.ind, kk)
                  }
              }
          }
      }
      ## Create segment test vector
      if(curr.nclust == prev.nclust + 1){
          curr.test.group.ind = (if(length(matched.curr.group.ind) == 0) (1:curr.nclust) else (1:curr.nclust)[-matched.curr.group.ind])
          group1 = curr.groups[[curr.test.group.ind[1]]]
          group2 = curr.groups[[curr.test.group.ind[2]]]
          newbasis[group1] = 1

          ## Calculate residual subspace basis for the added variable
          residual = resid(lm(newbasis ~ oldmodel-1))
          oldmodel = cbind(oldmodel, newbasis)

          models[[step+1]] = oldmodel
          residuals[[step+1]] = residual
      } else {
          models[[step+1]] = oldmodel
          residuals[[step+1]] = NA
      }

      ## update groups for next iteration
      prev.groups = curr.groups
      prev.nclust = curr.nclust
  }


  ## Calculate bic scores at dual knots
  bic = RSSs = complexities = rep(NA,maxsteps)
  for(ii in 1:length(df)){
    RSS = sum((resid(lm(y0 ~ models[[ii]]-1)))^2)
    complexity = log(length(y0))*(sigma^2)*df[ii]
    bic[ii] = RSS + complexity
    RSSs[ii] = RSS
    complexities[ii] = complexity
  }

  ## Also calculate bic at primal knots
  bic.primal = rep(NA,maxsteps+1)
  bic.primal[1] = bic[1]
  df.before = 1
  for(ii in 2:length(df)){
    if(df.before < df[ii]){
      RSS = sum((resid(lm(y0 ~ models[[ii]]-1)))^2)
#      complexity = log(df[ii])
      mypen = log(length(y0)) * (sigma^2) * df[ii]
      bic.primal[ii] = RSS + mypen#complexity
    }
    df.before = df[ii]
  }
  knots.primal = which(!is.na(bic.primal))
  models.primal    = models[!is.na(bic.primal)]
  residuals.primal = residuals[!is.na(bic.primal)]
  bic.primal       = bic.primal[!is.na(bic.primal)]

  return(list(models=models,df=df,residuals=residuals,bic=bic,
              bic.primal = bic.primal,
              knots.primal = knots.primal,
              models.primal,residuals.primal = models.primal,residuals.primal))
}


##' Helper function to get trend filtering contrast (for 'final model'
##' inference). Takes in |pathobject$action| and test time and stop time, and
##' returns the properly adjusted trendfiltering contrast
get.tf.contrast = function(action, y, test.time, stop.time, order=1){
    H = getH.trendfilter(n,order)
    test.coord = action[test.time] + 1
    adj.coord = c(1:(order+1), (1+states[[stop.time]]))
    adj.coord = adj.coord[adj.coord!=test.coord]
    X = H[, adj.coord]
    PH = X %*% solve(t(X) %*% X , t(X) )
    orthPH = (diag(n) - PH)
    v = as.numeric(H[,test.coord] %*% orthPH)
    if( v %*% y < 0) v = -v
    return(v)
  }



##' Helper for getdvec()
##' @param breaks all the breakpoints from path so far
##' @param k algorithm's step to use
##' @param klater later step to condition on
##' @param n length of response (= length of dvec )
makesegment  = function(breaks, k, klater, n){

if(length(breaks)<k) stop("not enough breaks!! k > number of breaks")
  K <- breaks[k] # is the index of the jump selected at step k

  # whether or not to condition on a later step (|klater|)
  kk <- klater

  relevantbreaks = (if(kk==1) c() else breaks[1:kk])
  endpoints = c(1,n)
  allbreaks <- c(endpoints, relevantbreaks)
  allbreaks <- sort(allbreaks)

  if(K %in% allbreaks) allbreaks = allbreaks[-which(allbreaks == K)]
  allbreaks = sort(unique(c(allbreaks, endpoints))) #temporary fix just in case the global endpoints are detected..
  min.index <- max(sum(allbreaks< K),1)             #temporary fix continued

  Kmin <- allbreaks[min.index]
  Kmax <- allbreaks[min.index + 1]

  if(Kmin != 1) Kmin = Kmin + 1 # special case handling

  return(list(Kmin=Kmin,K=K,Kmax=Kmax))
}




##' Takes in the signs |signs| and the locations |final.model| and returns a
##' vector of values that you can plot as a step sign plot.
##' @param signs signs of breakpoints
##' @param final.model breakpoint set
##' @param n length of data.
##' @examples
##' step.sign.plot(signs = c(+1,-1), final.model = c(20,36), n=60)
##' @export
step_sign_plot_inner = function(signs, final.model, n){
    signs = signs[order(final.model)]
    signs = c(-signs[1], signs)
    cumul.signs = cumsum(signs)
    final.model = sort(final.model)
    nn = length(final.model)
    indices <- vector(mode = "list", length = nn+1)
    if(length(final.model)==1){
        indices[[1]] = 1:final.model
        indices[[2]] = (final.model+1):n
    } else {
        indices[2:nn] = Map(function(a,b){a:b},
                            final.model[1:(length(final.model)-1)]+1,
                            final.model[2:length(final.model)])
        indices[[1]] = 1:final.model[1]
        indices[[nn+1]] = (final.model[length(final.model)]+1):n
    }

    sign0 = do.call(c,
                    lapply(1:length(cumul.signs),
                           function(ii){rep(cumul.signs[ii], length(indices[[ii]]))}))

    return(sign0)
}

