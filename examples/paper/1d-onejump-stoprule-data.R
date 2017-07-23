## Make sure you're working from [dropboxfolder]/code
source("settings.R")
source('funs.R')
source('examples/testfuns.R')
source('selectinf/selectiveInference/R/funs.inf.R')
source('dualPathSvd2.R')

## obtain things that will help calculate power
  n = 60 # 20,40,60,80 
  consec = 2
  nsim = 2000#50000
  ngrain = 15
  lev1 = 0;  lev2 = 1
  maxsteps = n-2
  sigma = 1
  lev2max = 4
  lev2list = seq(from=0, to = lev2max,length=ngrain)
  D = makeDmat(n,'trend.filtering',ord=0)
    

  # stoptimes & pval & verdict storage
  stoptimes.bic = stoptimes.ebic = stoptimes.aic = array(NA, c(ngrain,nsim))
  pvals.bic = pvals.bic.decluttered = pvals.ebic = pvals.aic = 
  verdicts.bic = verdicts.bic.decluttered = verdicts.ebic = verdicts.aic = array(NA,c(ngrain, nsim, n))
  pvals.fixed1 = pvals.fixed2 = verdicts.fixed1 = verdicts.fixed2 = array(NA,c(ngrain, nsim, n))
  pvals.oracle = verdicts.oracle = array(NA,c(ngrain,nsim))

  for(igrain in 1:ngrain){
    cat('\n',"one-jump simulations",'\n')
    cat('\n', igrain, "out of", ngrain, "for sample size" , n, "and maximum lev2", lev2max,'\n')
    lev2 = lev2list[igrain]
    for(isim in 1:nsim){
      cat('\r', isim, "out of", nsim)

      # Break if we have enough hits
      loc = n/2
      nhit.enough = 1000
      enough.hits.at.loc = function(verdict.array, igrain, loc, nhit.enough){ sum(!is.na(verdict.array[igrain,,loc])) > nhit.enough }
      if(all(unlist(lapply(list(verdicts.bic, verdicts.ebic, verdicts.aic), enough.hits.at.loc, igrain, loc, nhit.enough)))) {
        break
      }

      # Generate data + path
      beta0 = rep(c(lev1,lev2),each=n/2)
      y0    = beta0 + rnorm(n, 0,sigma)
      f0    = dualpathSvd2(y0,D,maxsteps,approx=T)
      
      mm = get.modelinfo(f0,y0,sigma,maxsteps,D=D,stoprule='bic')
      bic = mm$ic      
      stoptime.bic = which.rise(bic,consec) - 1 # internally defining the `stoptime' to be the step of the algorithm where you stop. the stoptime to be plotted is this+1.
      stoptime.bic = pmin(stoptime.bic, n-consec-1)
            
      if(stoptime.bic > 0){ 
      locs.bic = f0$pathobj$B[1:stoptime.bic]
      Gobj    = getGammat.with.stoprule(obj=f0,y=y0,condition.step=stoptime.bic+consec,
                                        type='tf',stoprule='bic',sigma=sigma,consec=consec,maxsteps=maxsteps, D=D)
      G       = Gobj$Gammat
      u       = Gobj$u
      
      # Non-decluttered
      for(test.step in 1:stoptime.bic){
        d       = getdvec(f0,y0,test.step,stoptime.bic,type="segment")
        loc     = locs.bic[test.step]
        pval    = pval.fl1d(y0, G, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE, u)
        verdict = (pval < (0.05/stoptime.bic))
        pvals.bic[igrain, isim, loc]       <- pval
        verdicts.bic[igrain, isim, loc]    <- verdict
      }
      
      # Decluttered
      states = get.states(f0$action)
      final.model = states[[stoptime.bic+1]]
      final.model.signs = f0$ss[[stoptime.bic+1]]
      final.model.decluttered = declutter(final.model)

      for(ii in 1:length(final.model.decluttered)){
        test.knot      = final.model.decluttered[ii]
        adj.knot       = final.model.decluttered
        test.knot.sign = final.model.signs[which(final.model == test.knot)]
        v =     make.v.tf.fp(test.knot = test.knot,
                             adj.knot  = adj.knot,
                             test.knot.sign = test.knot.sign,#f1$adj.knot,
                             D = D)        
        loc     = test.knot
        pval    = pval.fl1d(y0,G,dik=v,sigma,u=u)
        verdict = (pval < 0.05/length(final.model.decluttered))
        pvals.bic.decluttered[igrain, isim, loc]    = pval
        verdicts.bic.decluttered[igrain, isim, loc] = verdict
      }
      }


      # ebic: extended bic using factor of 1E-5
      ebic.fac=1E-3
      ebic = get.modelinfo(f0,y0,sigma,maxsteps,D=D,stoprule='ebic',ebic.fac=ebic.fac)$ic      
      stoptime.ebic = which.rise(ebic,consec) - 1 
      stoptime.ebic = pmin(stoptime.ebic,n-consec-1)
      if(stoptime.ebic > 0){
      locs.ebic= f0$pathobj$B[1:stoptime.ebic]
      Gobj    = getGammat.with.stoprule(obj=f0,y=y0,condition.step=stoptime.ebic+consec,
                                        type='tf',stoprule='ebic',sigma=sigma,consec=consec,maxsteps=maxsteps, D=D, ebic.fac=ebic.fac)
      G       = Gobj$Gammat
      u       = Gobj$u
      for(test.step in 1:stoptime.ebic){
        d = getdvec(f0,y0,test.step,stoptime.ebic,type="segment")
        loc = locs.ebic[test.step]
        pval    = pval.fl1d(y0, G, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE,u)
        verdict = (pval <  (0.05/stoptime.ebic))
        pvals.ebic[igrain, isim, loc]       <- pval
        verdicts.ebic[igrain, isim, loc]    <- verdict
      }
      }

      # aic
      aic = get.modelinfo(f0,y0,sigma,maxsteps,D=D,stoprule='aic')$ic 
      stoptime.aic = which.rise(aic,consec)-1
      stoptime.aic = pmin(stoptime.aic, n-consec-1)
      if(stoptime.aic > 0){ 
      locs.aic = f0$pathobj$B[1:stoptime.aic]
      Gobj    = getGammat.with.stoprule(obj=f0,y=y0,condition.step=stoptime.aic+consec,
                                      type='tf',stoprule='aic',sigma=sigma,consec=consec,maxsteps=maxsteps, D=D)
      G       = Gobj$Gammat
      u       = Gobj$u
      for(test.step in 1:stoptime.aic){
        d = getdvec(f0,y0,test.step,stoptime.aic,type="segment")
        loc     = locs.aic[test.step]
        pval    = pval.fl1d(y0, G,       d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE,u)
        verdict = (pval <  (0.05/stoptime.aic))
        pvals.aic[igrain, isim, loc]          <- pval
        verdicts.aic[igrain, isim, loc]       <- verdict
      }
      }
      
      # fixed stop times
      for( fixedstoptime in 1:2 ){
        G.truth = getGammat.naive(f0,y0,fixedstoptime)$G
        for( test.step in 1:fixedstoptime ){
          d = getdvec(f0,y0,test.step,fixedstoptime,type="segment")
          loc = f0$pathobj$B[test.step]
          pval = pval.fl1d(y0, G.truth, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE)
          verdict = (pval < (0.05/fixedstoptime))
          if ( fixedstoptime == 1 ){
            pvals.fixed1[igrain, isim, loc]    <- pval
            verdicts.fixed1[igrain, isim, loc] <- verdict
          } else if ( fixedstoptime == 2 ){
            pvals.fixed2[igrain, isim, loc]    <- pval
            verdicts.fixed2[igrain, isim, loc] <- verdict
          } else { print('not coded yet.. but this will never happen :)') }
        }
      }

      # store stoptimes
      stoptimes.bic[igrain,isim] = stoptime.bic
      stoptimes.ebic[igrain,isim] = stoptime.ebic
      stoptimes.aic[igrain,isim] = stoptime.aic
      
      # store oracle
      brk = n/2
      dif = abs(mean(y0[(brk+1):n]) - mean(y0[1:brk]))
      lvl = 0.05/2
      n1 = n2 = n/2
      z_crit = qnorm(1-lvl)*sigma*sqrt(1/n1 + 1/n2)
      verdicts.oracle[igrain,isim] = dif > z_crit
      pvals.oracle[igrain,isim]    = 1-pnorm(dif, mean=0, sd = sigma*sqrt(1/n1^2 + 1/n2^2))
    }
    
    # save
    obj.list1 = c("pvals.bic","pvals.bic.decluttered","pvals.ebic","pvals.aic",
                  "verdicts.bic","verdicts.bic.decluttered","verdicts.ebic","verdicts.aic",
                  "pvals.fixed1","pvals.fixed2","verdicts.fixed1","verdicts.fixed2",
                  "pvals.oracle","verdicts.oracle",
                  "stoptimes.bic","stoptimes.ebic","stoptimes.aic",
                   "sigma","nsim","ngrain","lev1","lev2","n", "lev2list","lev2max")#"sigmalist")
    save(list=obj.list1, file=file.path(outputdir, paste0("paper-BIC-onejump-segmentsize", n/2, ".Rdata")))
  }

