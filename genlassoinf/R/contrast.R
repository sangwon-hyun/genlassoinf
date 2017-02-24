### Function to make trend filtering segment contrasts
### from a regression point of view
#make.v.tf = function(test.knot, adj.knot, sign.test.knot, n, tf.order){
#  H = getH.trendfilter(n,1)
#  test.coord = test.knot+1
#  adj.coord = c(1:(tf.order+1), (1+adj.knot))
#  adj.coord = adj.coord[adj.coord!=test.coord]
#  X = H[, adj.coord ]
#  PH = fitted(lm(diag(n) ~ X -1))
#  orthPH = (diag(n) - PH)
#  v = as.numeric(H[,test.coord] %*% orthPH)
#  v = v * sign.test.knot
#  }
#  this.sign = f0$pathobj$s[f0$pathobj$B == final.model[ii]]
#  if(length(this.sign)==0){
#  f1 = dualpathSvd2(y0, D, maxsteps = stop.time-1)
#  this.sign = f1$pathobj$s[f1$pathobj$B == final.model[ii]]
#}




##' Function to make trend filtering segment contrasts
##' by use of null space projections
##' @param test.knot is the smaller model (knot set)
##' @param adj.knot is the larger model (knot set); \code{adj.knot} must be a superset of \code{test.knot}
##' @export
make.v.tf.fp = function(test.knot, adj.knot, test.knot.sign, D){# fp for first principles
    ## B1,B2 is smaller/larger model (knot sets)
    ## rr^T is pB1 - pB2, 
    ## pB1 formed by by proj null(D_{-B}) by orth proj onto row(D_{-B}) = col(t(D_{-B})) ~= tD

    stopifnot(test.knot %in% adj.knot)
    

  bigmodel = unique(c(adj.knot,test.knot)) #states[[ii]]
  smallmodel = c(adj.knot)
  smallmodel = smallmodel[smallmodel!=test.knot]
  tD = cbind(if(all(is.na(bigmodel))) t(D) else t(D)[,-(bigmodel)])
  proj = function(mymat){ return(mymat %*% solve(t(mymat)%*%mymat, t(mymat)))}

  # Big model projection
  rr = rankMatrix(tD)
  tDb = svd(tD)$u[,1:rr]
  big.proj = proj(tDb)

  # Small model projection
  tD.prev = cbind(if(all(is.na(smallmodel))) t(D) else t(D)[,-(smallmodel)])
  rr.prev = rankMatrix(tD.prev)
  tDb.prev = svd(tD.prev)$u[,1:rr.prev]
  small.proj = proj(tDb.prev)

  diff.proj = big.proj - small.proj
  myresid = svd(diff.proj)$u[,1]

  stopifnot(rankMatrix(diff.proj)==1 | as.numeric(rr.prev - rr) ==1)
  myresid = myresid / sqrt(sum((myresid)^2))
#  d1 = diff(myresid[(test.knot):(test.knot+1)])
#  d2 = diff(myresid[(test.knot+1):(test.knot+2)])
#  d2-d1
#  slopedir = d2-d1 > 0
#  slopedir = sign(slopedir-.5)

  slopedir = sign(D[test.knot,] %*% myresid)

  # This is ensuring that the direction of myresid at the test knot position (slopedir) 
  # is the same direction as the direction of slope change in the solution (test.knot.sign)
  # A conditional like # if(slopedir != sign.test.knot) .. # used to do the job
  
  myresid = test.knot.sign * slopedir  * myresid

  
  return(myresid)
}

##' In 1d fused lasso /only/, Get k'th contrast vector, for the segment or spike.
##' test.
##' @param obj either output from fl1d() or from dualpathSvd2()
##' @param y data vector (same as the one used to produce path)
##' @param k algorithm's step to use.
##' @param klater later step to condition on
get_v_1dfusedlasso = function(obj, y=NULL, k, klater = k, type =c("spike","segment"), n){

  type <- match.arg(type)
  if(k > klater) stop("Attempting to form contrast from a later step than the conditioning step.")
  
    ## ik and sk are the index and sign of the jump selected at step k
    ik = (obj$pathobjs)$B[k]
    sk = (obj$pathobjs)$s[k] 
    breaks = (obj$pathobjs)$B
    
    if(type == "spike"){
  
        v = rep(0,length(y))
        v[ik] = -1
        v[ik+1] = 1
        v = sk * v      
        
    } else if (type == "segment"){
        
        ## Extract usual segment test endpoints
        Ks = makesegment(breaks=breaks,k=k,klater=klater,n=length(y))
        K = Ks$K
        Kmin = Ks$Kmin
        Kmax = Ks$Kmax
        
        ## form vector
        v = rep(0,length(y))    
        v[Kmin:K] <- (Kmax - K)/(Kmax - Kmin + 1)
        v[(K+1):Kmax] <- -(K - Kmin + 1)/(Kmax - Kmin + 1)
        v <- -sk *v
  
    } else {

      stop("Not coded yet!")
        
    }
    return(v)
}

getdvec = get_v_1dfusedlasso

## Function to produce segment test contras ONLY for the fused lasso regression example in 2016 paper 
## TT is number of time points, J is number of stocks of interest
## type = "lrt" is for the likelihood ratio formulation, and type = "direct" is for forming segment test in 
## the coefficient scale first, then moving to the response scale
make.v.regression = function(test.knot, adj.knot, test.knot.sign, f0, TT, J,
                             X.orig,X.tilde, type=c("lrt","direct")){
    type = match.arg(type) 
    group.inds = lapply(1:J, function(ii){(TT*(ii-1)+1):(TT*(ii)-1) - (ii-1) })
    alt.breaks = unique(c(adj.knot , test.knot))#states[[stop.time]]
    null.breaks = alt.breaks[alt.breaks!=test.knot]    
  
    # Augment the design matrix at their breaks
    X.aug.alt = get.augmented.X(X.orig, alt.breaks,F,TT,J,group.inds) 
    X.aug.null = get.augmented.X(X.orig, null.breaks,F,TT,J,group.inds)
    
    ## # Ryan's contrast  
    ## v1 = rep(0,length(alt.breaks)+J)
    ## which.test.loc = which(test.knot == sort(alt.breaks))
    ## v1[which.test.loc] = -1 
    ## v1[which.test.loc+1] = 1
    ## coefs = solve(t(X.aug.alt)%*% X.aug.alt, t(X.aug.alt))
    ## v.direct = as.numeric(v1%*%(coefs))
    ## v.direct = v.direct * test.knot.sign
  
    # Ryan's contrast  
    v1 = rep(0,length(alt.breaks)+J)
    which.group = get.augmented.X(X.orig, alt.breaks,T,TT,J,group.inds)
    to.add = unlist(sapply(1:length(which.group), function(jj) rep(jj, length(which.group[[jj]]))))-1
    which.test.loc = which(test.knot == sort(alt.breaks))
    which.test.loc = which.test.loc + to.add[which.test.loc]
    print(paste("column in effective design matrix is", which.test.loc)) 
    v1[which.test.loc] = -1 
    v1[which.test.loc+1] = 1
    coefs = solve(t(X.aug.alt)%*% X.aug.alt, t(X.aug.alt))
    v.direct = as.numeric(v1%*%(coefs))
    v.direct = v.direct * test.knot.sign
    
    return(v.direct)
}


##' Different way of making regression segment test contrast 
make.v.regression.lrt = function(test.knot, adj.knot, test.knot.sign, f0, TT, J,
                             X.orig,X.tilde, type=c("lrt","direct")){
    type = match.arg(type) 
    group.inds = lapply(1:J, function(ii){(TT*(ii-1)+1):(TT*(ii)-1) - (ii-1) })
    alt.breaks = unique(c(adj.knot , test.knot))#states[[stop.time]]
    null.breaks = alt.breaks[alt.breaks!=test.knot]    
  
    # Augment the design matrix at their breaks
    X.aug.alt = get.augmented.X(X.orig, alt.breaks,F,TT,J,group.inds)
    X.aug.null = get.augmented.X(X.orig, null.breaks,F,TT,J,group.inds)
    
    # My contrast: the basis vector of the residual linear subspace
    P.alt = X.aug.alt %*% solve(t(X.aug.alt)%*% X.aug.alt, t(X.aug.alt))
    P.null = X.aug.null %*% solve(t(X.aug.null)%*% X.aug.null, t(X.aug.null))
    P.diff = P.null - P.alt
    if(rankMatrix(P.diff)!=1){
      print("rank of residual projection is not 1!")
      print(alt.breaks)
      print(null.breaks)
      print(rankMatrix(P.diff))
    } 
    my.contrast = svd(P.diff)$u[,1]
  
    ## Checking the signs using the immediately left/right adjacent segment means 
     translate = function(my.ind){ my.ind + floor(my.ind/(TT-1))}
     aug.breaks = sort(c(alt.breaks, (TT-1)*(0:J)))
     which.mid = which(test.knot == aug.breaks)
     which.left = which.mid - 1
     which.right = which.mid + 1
     right.inds = sapply((aug.breaks[which.mid]+1):aug.breaks[which.right], translate)
     left.inds = sapply((aug.breaks[which.left]+1):aug.breaks[which.mid], translate)
     left.mean = mean(my.contrast.coef.scale[left.inds])
     right.mean = mean(my.contrast.coef.scale[right.inds ])
     coef.scale.jump.direction = sign(right.mean-left.mean)
     my.contrast = coef.scale.jump.direction * test.knot.sign * my.contrast 
  
     return(my.contrast)
}



##' Helper function to /manually/ make contrasts.
##' @param test.bps Breakpoint location to test.
##' @param adj.bps Directly adjacent breakpoint locations.
##' @param sn Sign (direction) of proposed breakpoint at |test.bps|; use +1 or -1.
##' @param n length of data.
##' @examples
##' make_contrast(20,c(1,40),60)
##' @export
make_contrast = function(test.bp, adj.bps, sn, n){

    ## Basic checks
    stopifnot(all(c(test.bp, adj.bps) %in% 1:n))
    stopifnot(min(adj.bps)<=test.bp)
    stopifnot(max(adj.bps)>=test.bp)
    stopifnot(length(sn)==1)
    stopifnot(sn %in% c(-1,1))

    ## Make contrast
    d = rep(0,n)
    d[(min(adj.bps)):(test.bp)] = -1/(test.bp-min(adj.bps)+1)
    d[(test.bp+1):(max(adj.bps))] = +1/(max(adj.bps)-test.bp)
    return(sn*d)
}
