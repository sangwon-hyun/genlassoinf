##' From https://github.com/statsmaths/genlasso/blob/master/R/dualpathSvd.R
##' modified to compute the polyhedron post-selection inference selection.
##'
##' We compute a solution path of the generalized lasso dual problem:
##' \deqn{\hat{u}(\lambda) = \argmin_u \|y - D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty
##' \leq \lambda }
##'
##' where \eqn{D} is \eqn{m \times n}. Here there is no assumption on \eqn{D},
##' and we use a fresh SVD at each iteration (computationally naive but stable).
##'
##' Note: the df estimates at each \eqn{lambda_k} can be thought of as the df
##' for all solutions corresponding to lambda in \eqn{(lambda_k,lambda_{k-1})},
##' the open interval to the *right* of the current \eqn{lambda_k}.
##' @param y numeric vector of data.
##' @param D penalty matrix.
##' @param approx If approx=TRUE, then the fused lasso path will be computed
##'     without any coordinates leaving the path.
##' @param maxsteps maximum number of steps of the algorithms to run.
##' @param rtol Tolerance for solving linear system using svdsolve(), defaults
##'     to 1e-7.
##' @param ctol Tolerance for G&%y ~= c; defaults to 1e-7
##' @param btol Tolerance for leaving times, precision issue, defaults to 1e-10
##' @param cdtol Tolerance for cdtol; defaults to 1e-4
##' @param invalid temporary addition; getting rid of sign conditioning.
##' @return \code{path} class object, which is an R list containing results from
##'     having applied dual path algorithm. |nkstep| is a numeric vector that
##'     contains the row number of the polyhedron after each step.
##' @export
    
dualpathSvd2 <- function(y, D, approx=FALSE, maxsteps=2000, minlam=0,
                         rtol=1e-7, btol=1e-7, verbose=FALSE, object=NULL,
                         ctol=1e-10, cdtol=1e-4,
                         invalid = FALSE){

  # Error checking
  stopifnot(ncol(D) == length(y))

  nk = 0
  nkstep = 0
  shits.list = list()
  ss = list() ## list of all |s|'s along the way

  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)

    # Initialize Gamma matrix
    G = matrix(NA,nrow= maxsteps*ncol(D)*4 ,ncol=n)

    # Compute the dual solution at infinity, and find the first critical point.
    In = diag(1,n)
    sv = svdsolve(t(D),y,rtol)
    uhat = as.numeric(sv$x[,1])        # Dual solution
    q = sv$q                           # Rank of D

    ihit = which.max(abs(uhat))   # Hitting coordinate
    hit = abs(uhat[ihit])         # Critical lambda
    s = Sign(uhat[ihit])          # Sign
    k = 1

    if (verbose) {
      cat(sprintf("1. lambda=%.3f, adding coordinate %i, |B|=%i...",
                  hit,ihit,1))
    }

    # Now iteratively find the new dual solution, and
    # the next critical point

    # Things to keep track of, and return at the end
    buf = min(maxsteps,1000)
    u = matrix(0,m,buf)        # Dual solutions
    lams = numeric(buf)        # Critical lambdas
    h = logical(buf)           # Hit or leave?
    df = numeric(buf)          # Degrees of freedom
    action = numeric(buf)      # Action taken
    upol = c()                 # Constant in polyhedral constraint

    lams[1] = hit
    action[1] = ihit
    h[1] = TRUE
    df[1] = n-q
    u[,1] = uhat


    # add rows to Gamma
    tDinv = MASS::ginv(as.matrix(t(D)))

    # rows to add, for first hitting time (no need for sign-eligibility--just add all sign pairs)
    M = matrix(s*tDinv[ihit,], nrow(tDinv[-ihit,]), n, byrow=TRUE)
    newrows = rbind(M + tDinv[-ihit,],
                    M - tDinv[-ihit,])
    if(nrow(newrows)>=1) G[nk[length(nk)]+(1:nrow(newrows)),] = newrows
    nk = c(nk, nk[length(nk)]+nrow(newrows))
    nkstep[k] =  nk[length(nk)]

    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    Ds = D[ihit,]*s            # Vector t(D[B,])%*%s
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?


  } else {
  # If iterating an already started path
    # Grab variables needed to construct the path
    lambda = NULL
    for (j in 1:length(object)) {
      if (names(object)[j] != "pathobjs") {
        assign(names(object)[j], object[[j]])
      }
    }
    for (j in 1:length(object$pathobjs)) {
      assign(names(object$pathobjs)[j], object$pathobjs[[j]])
    }
    lams = lambda
  }

    ## Temporary!!!
    crits.list = list()
    crits = uhat
    names(crits) = Seq(1,m) ## not generalizable
    crits.list[[1]] = crits
    ## End of temporary

  tryCatch({
    while (k<=maxsteps && lams[k-1]>=minlam) {


      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        action = c(action,numeric(buf))
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        u = cbind(u,matrix(0,m,buf))
      }

      ##########
      # If the interior is empty, then nothing will hit
      if (r==m) {
        a = b = numeric(0)
        hit = 0
        q = 0
      }
      # Otherwise, find the next hitting time
      else {
        In = diag(1,n)
        sv = svdsolve(t(D1),cbind(y,Ds),rtol) #sv = svdsolve(t(D1),cbind(y,Ds,In),rtol)
        a = as.numeric(sv$x[,1])  # formerly a = as.numeric(D3 %*% y)
        b = as.numeric(sv$x[,2])
        D3 = MASS::ginv(as.matrix(t(D1)))# formerly as.matrix(sv$x[,3:(n+2)])
                                   # meant to be pseudoinverse of t(D[-I,])

        q = sv$q
        shits = Sign(a)
        shits.list[[k]] = shits ## temporary addition.
        hits = a/(b+shits);

        ## Temporary!!!
        crits = a/(b+shits)
        names(crits) = I
        crits.list[[k]] = crits
        ### End of temporary


        # Make sure none of the hitting times are larger
        # than the current lambda (precision issue)
        hits[hits>lams[k-1]+btol] = 0
        hits[hits>lams[k-1]] = lams[k-1]

        ihit = which.max(hits)
        hit = hits[ihit]
        shit = shits[ihit]

        # Gamma Matrix!
        # new rows to add, for viable hitting signs:
        tDinv = D3
        if(!invalid){
        newrows = (if(length(shits)>1){
          do.call(rbind, lapply(1:length(shits), function(ii){shits[ii] * tDinv[ii,]  }))
        } else {
          rbind(shits* tDinv)
        })
        } else {
        ## newrows = rbind(rep(NA,n))[-1,] #
        newrows = shits[ihit] * tDinv[ihit,]
        }

        ## Add those rows
        if(nk[length(nk)] + nrow(newrows) > nrow(G)){
            emptyrows = matrix(NA,nrow= nrow(G),ncol=n)
            G = rbind(G, emptyrows)
        }
        if(nrow(newrows)>=1) G[nk[length(nk)]+(1:nrow(newrows)),] = newrows
        nk = c(nk, nk[length(nk)]+nrow(newrows))

        # rows to add, for hitting event: (this is just dividing each row of
        # tDinv by corresponding element of tDinv%*%Ds+shit)
        A = asrowmat(D3/(b+shits)) # tDinv / as.numeric(tDinv %*% Ds + shits)
        if(nrow(A)!=1){

          nleft = nrow(A[-ihit,])
          if(is.null(nleft)) nleft = 1
          M = matrix(A[ihit,], nrow = nleft, ncol = n, byrow = TRUE)
          newrows = rbind(M - A[-ihit,])

          ## Add those rows
          if(nk[length(nk)] + nrow(newrows) > nrow(G)){
              emptyrows = matrix(NA,nrow= nrow(G),ncol=n)
              G = rbind(G, emptyrows)
          }
          if(nrow(newrows)>=1) G[nk[length(nk)]+(1:nrow(newrows)),] = newrows
            nk = c(nk, nk[length(nk)]+nrow(newrows))
        }

      }

      ##########
      # If nothing is on the boundary, then nothing will leave
      # Also, skip this if we are in "approx" mode
      if (r==0 || approx) {
        leave = 0
      }

      # Otherwise, find the next leaving time
      else {
        c = as.matrix(s*(D2%*%(y-t(D1)%*%a)))
        d = as.matrix(s*(D2%*%(Ds-t(D1)%*%b)))


        # round small values of c to zero (but not d)
        #cdtol = 1E-10
        c[abs(c) < cdtol] = 0

        # get leave times
        leaves = c/d

        # identify which on boundary set are c<0 and d<0
        Ci = (c < 0)
        Di = (d < 0)
        Ci[!B] = Di[!B] = FALSE
        CDi = (Ci & Di)

        # c and d must be negative at all coordinates to be considered
        leaves[c>=0|d>=0] = 0

        # Make sure none of the leaving times are larger
        # than the current lambda (precision issue)
        super.lambda = leaves>lams[k-1]+btol
        closeto.lambda = (lams[k-1] < leaves ) & (leaves  < lams[k-1]+btol)
        leaves[leaves>lams[k-1]+btol] = 0
        leaves[leaves>lams[k-1]] = lams[k-1]

        # If a variable just entered, then make sure it
        # cannot leave (added from lasso.R)
        if (action[k-1]>0) leaves[r] = 0

        # index of leaving coordinate
        ileave = which.max(leaves)
        leave = leaves[ileave]

        # Gamma matrix!!
        # rows to add, for leaving event:
        if(dim(D1)[1]==0) D1 = rbind(rep(0,ncol(D1))) # temporarily added because of
                                               # dimension problem in next line,
                                               # at last step of algorithm

        gmat = s*(D2%*%(In - t(D1)%*%D3)) # coefficient matrix to c
        # close-to-zero replacement is hard-coded in
        gmat[abs(c)<cdtol,] = rep(0,ncol(gmat))

        # we still want to see that gmat&%y ~= c
        if(!(max(gmat%*%y-c) < ctol)) print(max(gmat%*%y-c))

        gd = gmat / as.numeric(d)


        # hard coding in the zero replacements, to make it identical with lea
        gd[c>=0,] = rep(0, ncol(gd))        # re-doing zero-replacement in the g/d matrix
        gd[super.lambda,] = rep(0,ncol(gd)) # re-doing larger-than-lambda-replacement
        #gd[closeto.lambda,] = rep(0,ncol(gd)) # re-doing close-to-lambda-replacement (not sure how to make the i'th leaving time gd[i,]%*%y == lam[k-1] properly; solving to get some gd[i,] is like finding a non-unique solution to an overdetermined system; because such a gd[i,] is not unique, how do I know that adding this row to Gamma won't do weird and mysterious things?)

        #if( (length(Di)!=0) & (which(closeto.lambda) %in% which(Di))) print("closeto.lambda replacement SHOULD have happenned (but didn't).")

        # new rows that ensure c<0 #(only in )
        newrows = rbind(gmat[Ci&Di,]*(-1),
                         gmat[(!Ci)&Di,])

        ## Add those rows
        if(nk[length(nk)] + nrow(newrows) > nrow(G)){
            emptyrows = matrix(NA,nrow= nrow(G),ncol=n)
            G = rbind(G, emptyrows)
        }
        if(nrow(newrows)>=1) G[nk[length(nk)]+(1:nrow(newrows)),] = newrows
        nk = c(nk, nk[length(nk)]+nrow(newrows))

        # get rid of NA rows in G (temporary fix)
        missing.rows = apply(G, 1, function(row) any(is.na(row)))
        if(sum(missing.rows)>=1){ G <- G[-which(missing.rows),] }
        nk[length(nk)] = nrow(G)

        # new rows for maximizer
        CDi = (Ci & Di)
        CDi[ileave] = FALSE
        CDind = which(CDi)

        newrows = rbind(gd[rep(ileave,length(CDind)),] - gd[CDind,])

        ## Add those rows
        if(nk[length(nk)] + nrow(newrows) > nrow(G)){
            emptyrows = matrix(NA,nrow= nrow(G),ncol=n)
            G = rbind(G, emptyrows)
        }
        if(nrow(newrows)>=1) G[nk[length(nk)]+(1:nrow(newrows)),] = newrows
        nk = c(nk, nk[length(nk)]+nrow(newrows))


      }
      ##########
      # Stop if the next critical point is negative

      if (hit<=0 && leave<=0) {break}

      # If a hitting time comes next
      if (hit > leave) {

        # Record the critical lambda and solution
        lams[k] = hit
        action[k] = I[ihit]
        h[k] = TRUE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        u[,k] = uhat


        # add row to Gamma to characterize the hit coming next
        if(!approx)  {
            ## this is literally h_k - l_k > 0
            newrows = rbind(A[ihit,] - gd[ileave,])

            ## Add those rows
            if(nk[length(nk)] + nrow(newrows) > nrow(G)){
                emptyrows = matrix(NA,nrow= nrow(G),ncol=n)
                G = rbind(G, emptyrows)
            }
            if(nrow(newrows)>=1) G[nk[length(nk)]+(1:nrow(newrows)),] = newrows
            nk = c(nk, nk[length(nk)]+nrow(newrows))
        }

        # Update all of the variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        Ds = Ds + D1[ihit,]*shit
        s = c(s,shit)
        ss[[k]] = s
        D2 = rbind(D2,D1[ihit,])
        D1 = D1[-ihit,,drop=FALSE]
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      }

      # Otherwise a leaving time comes next
      else {
        # Record the critical lambda and solution
        lams[k] = leave
        action[k] = -B[ileave]
        h[k] = FALSE
        df[k] = n-q
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        u[,k] = uhat


        # new row to Gamma to characterize the leave coming next
        if(!approx)  {
          newrows = rbind(- A[ihit,] + gd[ileave,])

          ## Add those rows
          if(nk[length(nk)] + nrow(newrows) > nrow(G)){
              emptyrows = matrix(NA,nrow= nrow(G),ncol=n)
              G = rbind(G, emptyrows)
          }
          if(nrow(newrows)>=1) G[nk[length(nk)]+(1:nrow(newrows)),] = newrows
          nk = c(nk, nk[length(nk)]+nrow(newrows))
        }

        # Update all of the variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        Ds = Ds - D2[ileave,]*s[ileave]
        s = s[-ileave]
        ss[[k]] = s
        D1 = rbind(D1,D2[ileave,])
        D2 = D2[-ileave,,drop=FALSE]

        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }
      }

        nkstep[k] =  nk[length(nk)]

      # Increment step counter
      k = k+1
    }
  }, error = function(err) {
    err$message = paste(err$message,"\n(Path computation has been terminated;",
      " partial path is being returned.)",sep="")
    warning(err)})

    ## Trim
    lams = lams[Seq(1,k-1)]
    h = h[Seq(1,k-1)]
    df = df[Seq(1,k-1),drop=FALSE]
    u = u[,Seq(1,k-1),drop=FALSE]

    ## Save needed elements for continuing the path
    pathobjs = list(type="svd", r=r, B=B, I=I, approx=approx, k=k, df=df, D1=D1,
                    D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, h=h, rtol=rtol,
                    btol=btol, s=s, y=y)

    ## If we reached the maximum number of steps
    if (k>maxsteps) {
        if (verbose) {
            cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
            cat(" skipping the rest of the path.")
        }
        completepath = FALSE
    } else if (lams[k-1]<minlam) {
        ## If we reached the minimum lambda
        if (verbose) {
            cat(sprintf("\nReached the minimum lambda (%.3f),",minlam))
            cat(" skipping the rest of the path.")
        }
        completepath = FALSE
    } else {
        ## Otherwise, note that we completed the path
        completepath = TRUE
    }
    if (verbose) cat("\n")

    colnames(u) = as.character(round(lams,3))

    beta <-  apply(t(D)%*%u,2,function(column){y-column})
    states = get.states(action)
    ## cp = abs(action) ## This was wrong!!
    cp = states[[length(states)]]
    cp.sign = s#sign(action)

    obj = list(lambda=lams, beta=beta, fit=beta, hit=h, df=df, y=y, u=u, s=s, ss=ss,
               states=states, completepath=completepath,bls=y,
               pathobjs=pathobjs, nk = nk, nkstep = nkstep, action=action, D = D,
               cp = cp, cp.sign = cp.sign, ## Changepoint signs
               Gobj.naive=list(G=G[1:(nk[length(nk)]),],u=rep(0,nk[length(nk)])),
               Gobj.stoprule = NULL,
               maxsteps=maxsteps,
               crits.list=crits.list,## temporary!!!
               shits.list=shits.list)
    class(obj) = c("genlassoinf", "list")
    return(obj)
}


##' Function generic for stop_path()
##' @export
stop_path <- function(x,...) UseMethod("stop_path")


##' Embed stopping time into the \code{genlassoinf} class object \code{obj}, by
##' additionally inserting "stoptime", "stoprule", "consec", "stoppedmodel",
##' "Gobj.stoprule" into the object. Warning is issued if stoptime is zero, in
##' which case selected model is empty and the additional rows and entries to
##' \code{Gobj.stoprule$G} and \code{{Gobj.stoprule$u}} are empty as well.
##' @param obj A naive |path| object
##' @param sigma Standard deviation generating the data, in \code{obj$y}
##' @param stoprule Either one of \code{c("bic","ebic","aic")}. Defaults to
##'     \code{"bic"}.
##' @param consec How many rises do you want to stop at? Defaults to 2.
##' @examples
##' n = 60
##' sigma=1
##' consec = 2
##' D = makeDmat(n,type='tf',ord=0)
##' set.seed(1)
##' y0 = rep(c(0,1),each=n/2) + rnorm(n,0,sigma)
##' f0 = dualpathSvd2(y0,D,maxsteps,approx=T)
##' f0 = stop_path(f0)
##' @export
stop_path.genlassoinf <- function(obj, sigma, stoprule = "bic", consec = 2){

    ## Basic checks
    if(!is.null(obj$Gobj.stoprule)){
        stop("Stoprule already seems to have been applied. Original path object will be returned.")
        return(obj)
    }

    ## Get bic stopping time
    ic = get.modelinfo(obj, consec=2, sigma=sigma, stoprule = stoprule)$ic
    stoptime = which.rise(ic,consec) - 1
    stoptime = pmin(stoptime, length(obj$y)-consec-1)
    ## stopifnot(stoptime>0)
    if(stoptime<=0){
        warning("Stoptime is zero!")
        locs.bic = c()
        obj$Gobj.stoprule = list(G=rbind(rep(NA,length(obj$y)))[-1,], u=c())
        obj$stoppedmodel = obj$states[[stoptime+1]]
    } else {
        locs.bic = obj$pathobj$B[1:stoptime]

        ## Get object for stopping time.
        new.Gobj = getGammat.with.stoprule(obj=obj, y=obj$y,
                                       condition.step=stoptime+consec,
                                       type='tf', stoprule = stoprule, sigma=sigma,
                                       consec=consec, maxsteps=obj$maxsteps, D=D)
        obj$Gobj.stoprule = list(G=new.Gobj$G, u=new.Gobj$u)
        obj$stoppedmodel = obj$states[[stoptime+1]]
        }
    stoptime=0

    ## Add stopping time to path object
    obj$stoptime = matrix(stoptime,
                          dimnames=list(NULL,paste(stoprule, "with consec=", consec)))
    obj$stoprule = stoprule
    obj$consec = consec

    ## Return updated object
    class(obj) <- "genlassoinf"
    return(obj)
}


