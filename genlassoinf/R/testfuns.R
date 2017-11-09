##'  Function to run post-selection tests on graphs and store p-values results
##'  Currently only is able to do 3-clusters
##'  Input : path object, y0, initial graph, maximum steps to take, Dmat, cluster size, cluster #
##'  Output: 6 p values for each possible segment test.
sbmsim = function(f0,y0,mygraph, maxsteps,Dmat,clustersize,nclust=3,verbose=FALSE, stoptime=F){
  stopifnot(nclust==3)
  pvec = rep(NA,6)
  testtime.vec = rep(NA,6)

  # set initial graph's grouping
  prev.groups = list()
  for(member in unique(igraph::clusters(mygraph)$membership))  prev.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)

  prev.nclust = length(prev.groups)

  # do clustering and conduct tests along the way (and plot them)
  for(step in 1:maxsteps){


    if(verbose){  cat(step , "out of", maxsteps, fill = TRUE) }

    # Get information about that step
      G0 = getGammat(f0, y0, step)
      Dmat.curr = Dmat[-getB.from.actions(f0$action[1:step]),]
      mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat.curr), mode="undirected")
      if(nrow(rbind(Dmat.curr))==0) break


    # Keep track of connected components
      curr.groups = list()
      vsegment = rep(0,length(y0))

      for(member in unique(igraph::clusters(mygraph)$membership)){
        curr.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
      }

    curr.nclust = length(curr.groups)
    # if cluster was created
    if(curr.nclust == prev.nclust + 1){
      if(verbose==TRUE) { print("cluster created!") }
      matched.curr.group.ind = c()
      # Which cluster was created? Scan all groups from prev and curr; what has changed?
      for(jj in 1:prev.nclust){
        for(kk in 1:curr.nclust){
          if( setequal(prev.groups[[jj]], curr.groups[[kk]])){
            vsegment[prev.groups[[jj]]] = 0 # if _any_ of the previous groups match with current groups exactly, assign zero
            matched.curr.group.ind = c(matched.curr.group.ind, kk)
          }
        }
      }
      # create segment test vector
      curr.test.group.ind = (if(length(matched.curr.group.ind) == 0) (1:curr.nclust) else (1:curr.nclust)[-matched.curr.group.ind])
      group1 = curr.groups[[curr.test.group.ind[1]]]
      group2 = curr.groups[[curr.test.group.ind[2]]]
      vsegment[group1] = -1/length(group1)
      vsegment[group2] =  1/length(group2)


      # conduct segment test
      pseg = pval.fl1d(y0, G0, vsegment, sigma, approxtype="rob")

      # collect p-value only if it corresponds to one of the pairs!
      pair.list = list(list(1:clustersize, (1:clustersize)+clustersize ),
                      list((1:clustersize)+clustersize,(1:clustersize)+2*clustersize ),
                      list(1:clustersize,(1:clustersize)+2*clustersize),
                      list(1:(2*clustersize),(1:clustersize)+2*clustersize),
                      list(c(1:clustersize,(1:clustersize)+2*clustersize),(1:clustersize)+clustersize),
                      list(1:clustersize,(1+clustersize):(3*clustersize)))
      pair.list = lapply(pair.list,function(lst){ list(as.integer(lst[[1]]), as.integer(lst[[2]] ))})


      # record the p-value in the place where the test is defined!
      for(ll in 1:length(pair.list)){
         if( all(list(group1, group2) %in% pair.list[[ll]])){
           pvec[ll] = pseg
           testtime.vec[ll] = step
         }
      }

    } else if (curr.nclust +1 == prev.nclust){
      # if clusters was merged
       warning("cluster being CREATED has not been coded yet")
    } else {
      if(verbose==TRUE) { print("no cluster created!") }
    }

    # update groups for next iteration
    prev.groups = curr.groups
    prev.nclust = curr.nclust

    # break if there are enough clusters
#      if(precurr.nclust > 4){
#        cat("found enough clusters, exiting loop")
#        break
#      }
  }

  if(stoptime==T){
    return(list(pvec=pvec,testtime.vec=testtime.vec))
  } else {
    return(pvec)
  }
}





##' Runs post-selection tests on 2d graphs and store p-values , and outputs 6 p
##' values for each possible segment test.
##' @param f0 path object
##' @param y0 data
##' @param mygraph initial graph
##' @param maxsteps maximum steps to take
##' @param Dmat graph penalty matrix to be used
##' @param clustersize size of cluster
##' @param verbose whether to be loud or not
graph2dsim <- function(f0,y0,mygraph, maxsteps,Dmat,clustersize,verbose=FALSE){
    results = list()

    ## set initial graph's grouping
    prev.groups = list()
    for(member in unique(igraph::clusters(mygraph)$membership)){
        prev.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
    }

    prev.nclust = length(prev.groups)
    count = 1

    # Do clustering and conduct tests along the way (and plot them)
    for(step in 1:maxsteps){
      if(verbose){  cat(step , "out of", maxsteps, fill = TRUE) }

      # Get information about that step
      G0 = getGammat(f0, y0, step)
      Dmat.curr = Dmat[-getB.from.actions(f0$action[1:step]),]
      mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat.curr), mode="undirected")
      if(nrow(rbind(Dmat.curr))==0) break


      # Keep track of connected components
      curr.groups = list()
      vsegment = rep(0,length(y0))

      for(member in unique(igraph::clusters(mygraph)$membership)){
        curr.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
      }
      curr.nclust = length(curr.groups)

      # if cluster was created
      if(curr.nclust == prev.nclust + 1){
        if(verbose==TRUE) { print("cluster created!") }
        matched.curr.group.ind = c()
        # Which cluster was created? Scan all groups from prev and curr; what has changed?
        for(jj in 1:prev.nclust){
          for(kk in 1:curr.nclust){
            if( setequal(prev.groups[[jj]], curr.groups[[kk]])){
              vsegment[prev.groups[[jj]]] = 0 # if _any_ of the previous groups match with current groups exactly, assign zero
              matched.curr.group.ind = c(matched.curr.group.ind, kk)
            }
          }
        }
        # create segment test vector
        curr.test.group.ind = (if(length(matched.curr.group.ind) == 0) (1:curr.nclust) else (1:curr.nclust)[-matched.curr.group.ind])
        group1 = curr.groups[[curr.test.group.ind[1]]]
        group2 = curr.groups[[curr.test.group.ind[2]]]

        vsegment[group1] = +1/length(group1)
        vsegment[group2] =  -1/length(group2)

        # adjust the sign of vsegment
        if(vsegment%*%y0 <= 0) vsegment = -vsegment

        # conduct segment test
        pseg = pval.fl1d(y0, G0, vsegment, sigma, approxtype="rob")

        # record p-value and the groups that are being tested (and parse them later when aggregating results)
        results[[count]] = list(p = pseg, testttime = step, groups = list(group1=group1,
                                                            group2=group2))
        count = count + 1

      } else if (curr.nclust +1 == prev.nclust){
        # if clusters was merged
         warning("cluster being CREATED has not been coded yet")
      } else {
        if(verbose==TRUE) { print("no cluster created!") }
      }

      # update groups for next iteration
      prev.groups = curr.groups
      prev.nclust = curr.nclust
    }
    return(results)
}

##' Generates one jump data
onejump.y <- function(returnbeta=F, lev1=1, lev2=5, sigma=.5, n=60){
  beta = rep(c(lev1,lev2),each=n/2)
  if(returnbeta) return(beta)
  y <- beta + rnorm(n,sd=sigma)
  return(y)
}

twojump.y <- function(returnbeta=F,  n = 60,  sigma = .5){
  beta = rep(c(1,6,8),each=n/3)
  y <- beta + rnorm(n,sd=sigma)
  return(if(returnbeta) beta else y)
}

fourjump.y <- function(returnbeta = F,  n = 100,  sigma = 1,levs){
  if(length(levs)!=5) stop("Provide five levels!")
  beta = rep(levs,each=n/5)
  set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
  y = beta + rnorm(n,sd=sigma)
  return(if(returnbeta) beta else y)
}

alternjump.y <- function(returnbeta = F, lev1 = 1, lev2 = 5,  n = 60,  sigma = .5, seed = NA){
  beta = rep(c(lev1,lev2,lev1),each=n/3)
  if(is.na(seed)) {set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
  } else {set.seed(seed)}
  y = beta + rnorm(n,sd=sigma)
  return(if(returnbeta) beta else y)
}




##' Function to collect pvals, for Figure 8 and Appendix figure 8 and 15
onejump.naive.sim <- function(testtype = c("spike","segment"), nsim, sigma, lev1,
                              lev2list, numsteps, verbose=T, alpha=.05, loctype=c("exact", "oneoff"),
                              fac=30, gridsize=1000){

    testtype <- match.arg(testtype)
    loctype <- match.arg(loctype)

    if(loctype=="exact"){

        ## Initialize things
        pvals.correctlist = cis.correctlist = list()
        dlist = list()

        ## Collect p values when detection is correct
        for(ii in 1:length(lev2list)){
            lev2 = lev2list[ii]
            cat("\n", "lev2 is ", lev2, "\n")
            alldone = F
            jj = 0 ## Counter for how many hits we've had, on location lev/2
            pvals.correct = c()
            cis.correct = list()
            ds = list()
            while(!alldone){
                y = onejump(lev2,n) + rnorm(n,0,sigma)
                path   = dualpathSvd2(y,dual1d_Dmat(length(y)),maxsteps=numsteps,approx=T)
                G = path$Gobj.naive$G
                u = path$Gobj.naive$u
                d      = getdvec(obj=path, y=y, k=1, type=testtype, scaletype="segmentmean")
                if(path$pathobj$B[1] == n/2){
                    jj = jj+1
                    pvals.correct[jj] = poly.pval(y=y,G=G,v=d,u=u,sigma=sigma)$pv
                    cis.correct[[jj]] = confidence_interval(y, list(gamma=G,u=u), d,
                                                            sigma=1, alpha=alpha,
                                                            alternative="one.sided",
                                                            fac=fac, gridsize=gridsize)
                    ds[[jj]] = d
                    if(verbose) cat("\r", jj, "of", nsim)
                }
                alldone = (jj == nsim)
            }
            dlist[[ii]] = ds
            pvals.correctlist[[ii]] = pvals.correct
            cis.correctlist[[ii]] = cis.correct
        }

        return(list(pvals.correctlist = pvals.correctlist,
                    cis.correctlist = cis.correctlist,
                    dlist = dlist,
                    ## pvals.oneoff  = pvals.oneoff,
                    lev1 = lev1,
                    lev2list = lev2list,
                    nsim = nsim,
                    sigma = sigma,
                    alpha = alpha
                    ## kk = kk,
                    ## jj = jj,
                    ## prop.oneoff = kk/jj,
                    ))

    } else if (loctype == "oneoff") {

        ## collect p values from when detection is one off
        if(verbose) cat('\n','collecting p-values for 1 off','\n')
        alldone=F
        jj=ii=kk=0
        while(!alldone){
            y = onejump(lev2,n) + rnorm(n,0,sigma)
            path   = dualpathSvd2(y0,dual1d_Dmat(n),maxsteps=numsteps,approx=T)
            G = path$Gobj.naive$G
            u = path$Gobj.naive$u
            d      = getdvec(obj=path, y=y0, k=1, type=testtype, scaletype="segmentmean")
            if(abs(path$pathobj$B[1] - n/2) == 1){
                pvals.oneoff[jj] = poly.pval(y=y0,G=G,v=d,u=u,sigma=sigma)$pv
                jj = jj+1
                if(verbose) cat("\r", jj, "of", nsim)
            }
            kk = kk+1
            alldone = (jj == nsim)
        }

        return(list(pvals.oneoff = pvals.oneoff,
                    lev1 = lev1,
                    lev2list = lev2list,
                    nsim = nsim,
                    sigma = sigma,
                    alpha = alpha,
                    kk = kk,
                    jj = jj,
                    prop.oneoff = kk/jj
                    ))

    } else {
        stop(paste(loctype, "type of |loctype| isn't handled!"))
    }
}


onejump <- function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}





##' Function to aggregate successes and hits in order to calculate condit
##' powers, for stopping time simulations
getpowers.from.chunks <- function(n, verdict.obj.name, ngrain= 20,
                                  nsim = 100000, nchunks = 100,
                                  file.id = "bic-onejump-segmentsize-allbics"){

    # Helper functions
    getpow = function(verdicts,ii,loc){  return(sum(verdicts[ii,,loc],na.rm=T)/pmax(1,sum(!is.na(verdicts[ii,,loc]))))   }
    getpow.numer = function(verdicts){  return( sum(verdicts,na.rm=T))   }
    getpow.denom = function(verdicts){  return( pmax(1,sum(!is.na(verdicts))))   }

    # Powers
    powers = array(NA, c(ngrain,n))
    powers.numer = powers.denom = array(0,c(ngrain,n))
    powers.prox.numer = powers.prox.denom = powers.prox = rep(0,ngrain)

    loc = n/2
    proxlocs = (loc-log(n)):(loc+log(n))

    # extract numer and denom and sum over chunks
    for(chunk.i in 1:nchunks){
    cat('\r', chunk.i, "out of", nchunks)
    load(file=file.path(codedir, paste0("maxoutput/",file.id , n/2, "-chunk",chunk.i,".Rdata")))
#      verdict.obj.name = "verdicts.bic2"
      verdicts = eval(parse(text=verdict.obj.name))
      for(ii in 1:ngrain){
        # Exact powers
        for(loc in 1:n){
          powers.numer[ii,loc]  = powers.numer[ii,loc] + getpow.numer(verdicts[ii,,loc])
          powers.denom[ii,loc]  = powers.denom[ii,loc] + getpow.denom(verdicts[ii,,loc])
        }
        # Approximate powers at true break coordinate
        prox.verdicts = apply(verdicts[ii,,proxlocs], 1,
                            function(this.sim.verdicts){
                             (if(!all(is.na(this.sim.verdicts))) any(this.sim.verdicts,na.rm=T) else NA)})
        powers.prox.numer[ii]  = powers.prox.numer[ii] + getpow.numer(prox.verdicts)
        powers.prox.denom[ii]  = powers.prox.denom[ii] + getpow.denom(prox.verdicts)
      }

      powers.prox.numer.bic = powers.prox.numer
      powers.prox.denom.bic = powers.prox.denom
    }
    cat('\n')

    # calculate condit power at each locaiton from extracted/summed numer and denom
    for(ii in 1:ngrain){
      powers[ii,]  = powers.numer[ii,]/powers.denom[ii,]
      powers.prox[ii] = powers.prox.numer[ii]/powers.prox.denom[ii]
    }
    rownames(powers) = names(powers.prox) = sigmalist
    colnames(powers) = 1:n

    return(list(powers=powers, powers.numer=powers.numer, powers.denom=powers.denom, powers.prox = powers.prox))
}
