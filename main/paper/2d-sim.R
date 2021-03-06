## Make sure you're working from [dropboxfolder]/code
  source("settings.R")
  source('selectinf/selectiveInference/R/funs.inf.R')
  source('funs.R')
  source('testfuns.R')
  source('dualPathSvd2.R')
  library(genlasso)
  library(igraph)
  verbose = F

## Generate simulation results
    nsim = 100000
    nhit.enough = 3000
    lev1 = 0
    lev2list = seq(from=15,to=75,by=15)/15 #  seq(from=15,to=75,by=15)/
    lev2list = c(lev2list, 0)
    ngrain = 5+1
    sigma=15/15
    nn = 10
    Result = list()
    Dmat = graph2D_Dmat(nn^2)
    maxsteps = 40
    consec = 1
    mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat), mode="undirected")
    chunks = list(1,2,3,4:5,6)
    nchunk = 4
    ## ichunk = 1
    latecounter = 0
    latestoptimes = rep(NA,nsim)

  pvals = stoptimes = primal.numchanges = matrix(NA,nrow=nsim,ncol=ngrain)
    for(igrain in chunks[[ichunk]]){#1:ngrain){
        lev2 = lev2list[igrain]
        cat(igrain, "out of", ngrain, "with nn x nn image, nn= ", nn, "for lev2=",lev2, fill=T)
        maxsteps = 50-igrain*4 #min(40+igrain*10,70)
        if(igrain==6) maxsteps=100 
        how.many.hits = 0
        for(isim in 1:nsim){
        cat("\r", isim, "out of", nsim, "with", how.many.hits, "hits so far, aiming for", nhit.enough,"hits")
        if(how.many.hits>nhit.enough) break

        # Make a noisy nn x nn image
            A = matrix(lev1,ncol=nn, nrow = nn)
            A[1:(nn/2),1:(nn/2)] = lev2
            A.noisy = A + rnorm(nn^2,mean=0,sd=sigma)

        # generate data
            beta0 = do.call(c,lapply(1:nrow(A),function(irow)A[irow,]))
            y0 = do.call(c,lapply(1:nrow(A.noisy),function(irow)A.noisy[irow,]))
            # optionally, plot!
            image(A.noisy)

        # fit graph fused lasso using the matrix given by D only (initially go a generous number of steps)
            f0 = dualpathSvd2(y0,Dmat,verbose=FALSE, maxsteps = maxsteps)
            maxsteps = length(f0$action[f0$action!=0])

        # Get BIC scores
            mm = get.modelinfo(f0,y0,sigma,maxsteps,Dmat, stoprule = 'bic')
            bic = mm$ic
            stop.time = which.rise(bic,consec)-1
            if(stop.time >= maxsteps){#print(paste("didn't stop by step", maxsteps));
                latecounter = latecounter+1
                latestoptimes[ii] = stop.time
                next
            }

        # Get polyhedron (with BIC)
            if(stop.time+consec < maxsteps){
            Gobj.new.with.stoptime = getGammat.with.stoprule(obj=f0,y=y0,
                                            condition.step = stop.time+consec,
                                            stoprule = "bic", sigma=sigma,
                                            consec=consec, maxsteps=maxsteps, type = "graph", D = Dmat)
            } else {
            print('bic rule hasnt stopped!')
            next
            }
            G = Gobj.new.with.stoptime$Gammat
            u = Gobj.new.with.stoptime$u

            if(!polyhedron.checks.out(y0,G,u)){
                print("polyhedron is problematic")
                next
            }

        # Calculate p-value and store it
        steps.to.test = mm$knots.primal[mm$knots.primal<=stop.time]
        steps.to.test = steps.to.test[steps.to.test!=1]
        states = get.states(f0$action)
        final.model = states[[stop.time+1]]
        final.model.signs = f0$ss[[stop.time+1]]

        for( mystep  in steps.to.test ){

            # Get the residual subspace basis
            v = mm$resids[,mystep]

            # Adjust it by the sign
            test.knot = final.model[mystep-1]
            slopedir = sign(Dmat[test.knot,] %*% v)
            test.knot.sign = final.model.signs[mystep-1]
            v = v * test.knot.sign * sign(slopedir-.5) 

            if(any(is.na(v))) print(mm)

            contrast.is.right=function(v){
            factors = factor(beta0)
            unique.factors = levels(factor(beta0))
            ind.1 = which(factors == unique.factors[1])
            ind.2 = which(factors == unique.factors[2])
            v = unname(v)
            return( (all.equal(v[ind.1],rep(v[ind.1[1]], length(ind.1)))==T) &
                    (all.equal(v[ind.2],rep(v[ind.2[1]], length(ind.2)))==T) )
            }

            if(contrast.is.right(v)){         
            how.many.hits = how.many.hits + 1
            pvals[isim,igrain] = pval.fl1d(y0, G, v, sigma, u = u)
            stoptimes[isim,igrain] = stop.time
            primal.numchanges[isim,igrain] = length(steps.to.test)

            if(is.nan(pvals[isim,igrain])) stop("nan pvalue occurred")
            }
        }
        }
        cat(fill=T)
        save(ngrain,nsim,lev2list,Dmat,maxsteps,mygraph,pvals, stoptimes, igrain, isim, primal.numchanges,chunks,ichunk,
            file = file.path(outputdir, paste0("paper-2dgraph-chunk-",ichunk,".Rdata")))
    }
cat(fill=T)
cat(latecounter, "runs out of", nsim, "stopped after |maxsteps|=", maxsteps, fill=T)
cat("and the stop times were, on average", mean(latestoptimes), "with std error", sd(latestoptimes),fill=T)
## cat("in fact, let me show them all to you! Here they are:", latestoptimes)
