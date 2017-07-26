##' Gets k'th Gammat matrix from path object \code{obj}. Note! If you're being
##' naive, then do whatever you want with condition.step.  Note, if you're using
##' stop times, condition.step should be the (stoptime + 2)
##' @export
getGammat.naive = function(obj, y, condition.step=NULL){

    ## Error checking
    n = length(y)
    if(length(obj$action) < condition.step ){
        stop("\n You must ask for polyhedron from less than ",
             length(obj$action), " steps! \n")
    }

    ## Extract naive Gamma matrix
    my.Gammat = obj$Gobj.naive$G[1:obj$nk[condition.step+1],] ## nk starts with a 0 element.
    my.u = rep(0,nrow(my.Gammat))
    return(list(G = my.Gammat, u = my.u, condition.step=condition.step))
}

##' Produces the rows to add to Gammat matrix for IC-based stopping rule in
##' generalized lasso signal approximation problem.
##' @export
getGammat.stoprule = function(obj, y, condition.step, stoprule, sigma, consec,
                              maxsteps, ebic.fac= 0.5, D, verbose=F){

    ## Basic checks
    if(!(stoprule %in% c('bic','ebic','aic'))){    stop('stoprule not coded yet!') }

    n = length(y)
    mm = get.modelinfo(obj=obj, consec=2, sigma=sigma, stoprule=stoprule,
                       ebic.fac=ebic.fac, verbose=verbose)

    ## get IC vector
    ic = mm$ic #get.ic(y,obj, sigma,maxsteps,D,type=stoprule,ebic.fac=ebic.fac)$ic

    ## s* : change in model size
    actiondirs = c(NA, sign(obj$action)[1:(maxsteps-1)])
    ## s** : change in bic
    seqdirs = c(getorder(ic))
    ## multiply s* and s**
    signs = seqdirs*actiondirs


    ## get stop times
    stoptime = which.rise(ic, consec, "forward") - 1

    if(length(seqdirs) < stoptime + consec + 1) stop(paste(stoprule, "has not stopped"))

    ## Check if conditioning steps are the same (i.e. stoptime + consec)
    stopifnot(condition.step == stoptime+consec)

    ## safely making enough empty rows to add
    rows.to.add = more.rows.to.add = matrix(NA, nrow = (stoptime+consec),
                                            ncol = ncol(obj$Gobj.naive$G))
    upol = more.upol = rep(NA,stoptime+consec-1)
    states = get.states(obj$action)

    ## Add rows to |G| and |u|
    for(jj in 1:(stoptime+consec)){

      residual = mm$resids[,jj+1]
      const    = abs(mm$pen[jj+1] - mm$pen[jj]) #getconst(stoprule,jj,sigma,n,big.df,small.df,bic.fac,ebic.fac)
      upol[jj] = (-signs[jj+1]) * sqrt(const)
      rows.to.add[jj,] = (-signs[jj+1]) * sign(t(residual)%*%y) * residual/sqrt(sum(residual^2))

      # Also add other direction if (1) model size change and (2) BIC comparison
  #    # are *both* up or down, in direction
  #    if(actiondirs[jj+1] == seqdirs[jj+1] & !any(is.na(residual))){
  #      more.upol[jj] = (-signs[jj+1]) * sqrt(const)
  #      more.rows.to.add[jj,] = (signs[jj+1]) * sign(t(residual)%*%y) * residual/sqrt(sum(residual^2))
  #    }
    }

    rows.to.add = rbind(rows.to.add,
                        more.rows.to.add)
    upol = c(upol, more.upol)

    # Get rid of non-primal-changing comparisons
    upol = upol[which(!is.na(rows.to.add[,1]))]
    rows.to.add = rows.to.add[which(!is.na(rows.to.add[,1])),]

    return(list(G = rows.to.add, u = upol))
}


##' Produces the rows to add to Gammat matrix, in the regression case requires
##' X.orig, ginvX.orig and D.orig which are from the pre-elastic-net regression
##' problem
##' @export
getGammat.stoprule.regression = function(obj, y, condition.step, stoprule,
                                           sigma, consec, maxsteps,  X.orig,
                                           ginvX.orig=NULL, D.orig,y0.orig){
    if(stoprule!='bic'){
      stop("Regression stopping rule only coded for BIC!")
    }

    n = length(y)
    if(is.null(ginvX.orig)) ginvX.orig = ginv(X.orig)

    # Function to get the entire sequence of information criterion values
    ic = getbic.regression(y0.orig,obj,sigma,maxsteps=maxsteps-1,
                X.orig = X.orig, ginvX.orig = ginvX.orig, D.orig = D.orig)


    # s* : change in model size
    actiondirs = c(NA, sign(obj$action)[1:(maxsteps-1)])
    # s** : change in bic
    seqdirs = c(getorder(ic))


    # multiply s* and s**
    signs = seqdirs*actiondirs

    # get stop times
    stoptime = which.rise(ic, consec, "forward") - 1

    if(length(seqdirs) < stoptime + consec + 1) stop(paste(stoprule, "has not stopped"))

    # Check if conditioning steps are the same (i.e. stoptime + consec)
    stopifnot(condition.step == stoptime+consec)

    # safely making enough empty rows to add
    rows.to.add = matrix(NA, nrow = (stoptime+consec),
                         ncol = ncol(obj$Gobj.naive$G))
    upol = rep(NaN,stoptime+consec-1)

    states = get.states(obj$action)

    # Add rows to Gamma matrix and |u|)
    # First comparison is the one between the null model and the 1st dual addition
  #  end.step = pmin( stoptime+consec, maxsteps)
    group.inds = lapply(1:J, function(ii){(TT*(ii-1)+1):(TT*(ii)-1) - (ii-1) })
    for(jj in 1:(stoptime+consec)){

      # Identify small and big model from which to harvest the residual
      prev.state = states[[jj]]
      this.state = states[[jj+1]]

      big.model = (if(actiondirs[jj+1]==+1){this.state} else {prev.state})
      small.model = (if(actiondirs[jj+1]==+1){prev.state} else {this.state})
      this.hit = big.model[!(big.model %in% small.model)]

      this.hit = obj$action[jj+1]
      alt.breaks = big.model#unique(c(states[[jj]] , this.hit))
      null.breaks = small.model#alt.breaks[alt.breaks!=this.hit]

      # Function to obtain the "augmented" X matrix, by breaking at the fused lasso breakpoints
      get.augmented.X = function(X,breaks,return.groups = F){
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

      # Augment the design matrix at their breaks
      X.aug.alt  = get.augmented.X(X.orig, alt.breaks )
      X.aug.null = get.augmented.X(X.orig, null.breaks)


      # projection function
      projection = function(mat){
        rm = invisible(rankMatrix(mat))
        if(rm<ncol(mat)){
          b = svd(mat)$u[,1:rm]
          return(b %*% solve(t(b)%*%b, t(b)))
        } else {
        return(mat %*% solve(t(mat)%*%mat, t(mat)))
        }
      }

      # Get the basis vector of the residual linear subspace
      P.alt = projection(X.aug.alt)#X.aug.alt %*% solve(t(X.aug.alt)%*% X.aug.alt, t(X.aug.alt))
      P.null = projection(X.aug.null)#X.aug.null %*% solve(t(X.aug.null)%*% X.aug.null, t(X.aug.null))
      P.diff = P.null - P.alt
      if(invisible(rankMatrix(P.diff))!=1){
        print("rank of residual projection is not 1!")
        print(alt.breaks); print(null.breaks);  print(rankMatrix(P.diff))
      }
      residual = svd(P.diff)$u[,1]

      # Get the constant
      const = (sigma^2)*(log(n))

      # Add it
      upol[jj] = (-signs[jj+1]) * sqrt(const)
      rows.to.add[jj,] = (-signs[jj+1]) * sign(t(residual)%*%y0.orig) * residual/sqrt(sum(residual^2))
    }

    # Get rid of non-primal-changing comparisons
    upol = upol[which(!is.na(rows.to.add[,1]))]
    rows.to.add = rows.to.add[which(!is.na(rows.to.add[,1])),]

    return(list(G = rows.to.add, u = upol))
}



##' Function that takes the same arguments as the getGammat.naive(), and some
##' extra arguments to extract the conditioning for the stopping rule.  Returns
##' the final Gamma matrix and u vector with stopping rules incorporated.
##' @export
getGammat.with.stoprule = function(obj, y, condition.step,# usage = c("fl1d", "dualpathSvd"),
                                   stoprule, sigma, maxsteps, consec,
                                   ebic.fac= 0.5,
                                   type = c("tf","graph", "regression"), X=NULL,
                                   D=NULL, X.orig=NULL, ginvX = NULL,
                                   D.orig=NULL, y0.orig = NULL,
                                   ginvX.orig = NULL,verbose=F){

    ## Basic error checking
    type = match.arg(type)
    if(type=="regression" & (is.null(y0.orig)|is.null(X.orig)|is.null(D.orig))){
        stop("You need to supply y, design matrix and D matrix from pre-elastic-net regression problem!")
    }

    ## Get naive stuff
    G.naive = getGammat.naive(obj, y, condition.step)#, usage = c("fl1d", "dualpathSvd"))

    ## Get stoprule-related stuff
    if(type %in% c("tf", "graph")){
        G.new.obj = getGammat.stoprule(obj=obj, y=y,
                                       condition.step=condition.step,
                                       stoprule=stoprule, sigma=sigma,
                                       consec=consec, maxsteps=maxsteps,
                                       ebic.fac= ebic.fac, D=D, verbose=verbose)
    } else if (type == "regression"){
        G.new.obj = getGammat.stoprule.regression(obj=obj,y=y,
                                                  y0.orig = y0.orig,
                                                  condition.step = condition.step,
                                                  stoprule = "bic", sigma=sigma,
                                                  consec=consec,
                                                  maxsteps=maxsteps,
                                                  X.orig = X.orig,
                                                  D.orig = D.orig,
                                                  ginvX.orig=ginvX.orig)
    } else {
        stop(paste(type, "not coded yet."))
    }

    ## Combine naive rows and new rows
    G.combined = rbind(G.naive$G,
                       G.new.obj$G)
    uvec = c(rep(0,nrow(G.naive$G)),
             G.new.obj$u)

    return(list(G = G.combined, u = uvec))
}



