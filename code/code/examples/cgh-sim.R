# working directory should be [genlassoinf]/code
#workingdir = '~/../../media/shyun/Bridge/Dropbox/CMU/courses(CURRENT)/genlassoinf/code'
#setwd(workingdir)
source("settings.R")
source('funs.R')
source('examples/testfuns.R')
source('dualPathSvd2.R')
lapply(c("genlasso","pryr"), require, character.only = TRUE)
source('selectinf/selectiveInference/R/funs.inf.R')
#library(polypath)

# call data
  library(cghFLasso)
  data(CGH)
  summary(CGH)
  
###########################################
######## Sparse Fused Lasso ###############
###########################################

# Fit sparse fused lasso paths
  y0 = CGH$GBM.y
  n = length(y0)
  maxsteps= 35#100
  D = rbind(makeDmat(n,order=0),diag(rep(1,n)))
  f0 = dualpathSvd2(y0,D,maxsteps=maxsteps,approx=T,verbose=T)
  tf.order = 0; consec = 1
  sigma = 0.46 # tuned in paper-signpic-cgh.R
  modelinfo = get.modelinfo(f0,y0,sigma,maxsteps,D,'bic',verbose=T)
  bic = modelinfo$ic#get.modelinfo(f0,y0,sigma,maxsteps,D,'bic',verbose=T)$ic

# Setup for model fitting
# Get stopping time and Polyhedron
  consec=1
  stop.time = which.rise(bic, consec=consec) - 1
  plot(bic); abline(v=stop.time+1, lty=2)
  Gobj = getGammat.with.stoprule(obj=f0,y=y0,
                                     condition.step = stop.time + consec,
                                     stoprule = "bic", sigma=sigma, ebic.fac=ebic.fac, type = "tf",
                                     consec=consec, maxsteps=maxsteps,D=D,verbose=T)
#  save(y0, f0, maxsteps, n, modelinfo, bic, D, consec, sigma, G, u, stop.time, file=file.path(outputdir,"cgh.Rdata"))

  load(file=file.path(outputdir,"cgh.Rdata")) # this is the pure-sparsity-inducing version
  G = Gobj$Gammat
  u = Gobj$u
   

# Conduct tests and record pvals + stoptimes
  states = get.states(f0$action)
        
# Declutter the last states
  final.model = states[[stop.time+1]]
  #which.primalchange = final.model<length(y0)
  #final.model = final.model[which.primalchange]
  final.model.signs = f0$ss[[stop.time+1]]   #[which.primalchange]
  ord = order(final.model)

  final.model.cluttered = final.model[ord]
  final.model.signs = final.model.signs[ord]
  final.model.early = final.model[final.model <= n]
  final.model.late  = final.model[final.model >  n]
  nbhood=2
  final.model.decluttered = c(declutter(final.model.early,nbhood), final.model.late)
  final.models = list(final.model.cluttered,final.model.decluttered)
#  final.models =  lapply(final.models,function(mymodel) mymodel[mymodel<=n])

  pvals = pvals.decluttered = array(NA,dim=c(n,maxsteps))
  results = results.decluttered = list()
  for(kk in 1:2){
    cat("\n type",kk,fill=T)
    final.model = final.models[[kk]]
    if(stop.time > 0){
      ntests = length(final.model)
      for(ii in 1:ntests){
        test.knot      = final.model[ii]        
        adj.knot       = final.model
        test.knot.sign = final.model.signs[which(final.model.cluttered == test.knot)]
        if(test.knot > n) next()
        print(test.knot)  

        
        # Check if there has been a change in primal
        bigmodel = unique(c(adj.knot,test.knot)) #states[[ii]]
        smallmodel = c(adj.knot)
        smallmodel = smallmodel[smallmodel!=test.knot]
        tD = cbind(if(all(is.na(bigmodel))) t(D) else t(D)[,-(bigmodel)])
        tD.prev = cbind(if(all(is.na(smallmodel))) t(D) else t(D)[,-(smallmodel)])
        if( rankMatrix(tD) - rankMatrix(tD.prev) == 0 ) {print("no change in primal, moving on"); next()}

        # make the differences based on the _reduced_ D matrix
        D.reduced = makeDmat(n,order=0)
        adj.knot = adj.knot[adj.knot<=n]
        print(adj.knot)

        # Calculate contrast vector
        v =     make.v.tf.fp(test.knot = test.knot,
                            adj.knot  = adj.knot,
                            test.knot.sign = test.knot.sign,#f1$adj.knot,
                            D = D.reduced)#D=D
        
        coord = final.model[ii]
        pval  = pval.fl1d(y0,G,dik=v,sigma,u=u)
        print(pval)

        if(kk==1){
          results[[ii]] = list(pval=pval,coord=coord,v=v)
        } else {
          results.decluttered[[ii]] = list(pval=pval,coord=coord,v=v)
        }
      }
    }    
  }
  save(file=file.path(outputdir,"cgh-results.Rdata"), 
       y0, f0, maxsteps, n, modelinfo, bic, D, consec, sigma, G, u, stop.time,results, results.decluttered, final.model.cluttered, final.model.decluttered)



###############################################
######## Non Sparse Fused Lasso ###############
###############################################

# Fit sparse fused lasso paths
  y0 = CGH$GBM.y
  n = length(y0)
  maxsteps= 100
  D = makeDmat(n,order=0)
  f0 = dualpathSvd2(y0,D,maxsteps=maxsteps,approx=T,verbose=T)
  tf.order = 0; consec = 1
  sigma = 0.46 # tuned in paper-signpic-cgh.R
  modelinfo = get.modelinfo(f0,y0,sigma,maxsteps,D,'bic',verbose=T)
  bic = modelinfo$ic#get.modelinfo(f0,y0,sigma,maxsteps,D,'bic',verbose=T)$ic

# Get stopping time and Polyhedron
  consec=2
  stop.time = which.rise(bic, consec=consec) - 1
#  plot(bic); abline(v=stop.time+1, lty=2)
  Gobj = getGammat.with.stoprule(obj=f0,y=y0,
                                     condition.step = stop.time + consec,
                                     stoprule = "bic", sigma=sigma, ebic.fac=ebic.fac, type = "tf",
                                     consec=consec, maxsteps=maxsteps,D=D,verbose=T)
  G = Gobj$Gammat
  u = Gobj$u
  save(y0, f0, maxsteps, n, modelinfo, bic, D, consec, sigma, G, u, stop.time, file=file.path(outputdir,"cgh-nonsparse.Rdata"))

  load(file=file.path(outputdir,"cgh-nonsparse.Rdata"))

# Conduct tests and record pvals + stoptimes
  states = get.states(f0$action)
        
# Declutter the last states
  final.model = states[[stop.time+1]]
  #which.primalchange = final.model<length(y0)
  #final.model = final.model[which.primalchange]
  final.model.signs = f0$ss[[stop.time+1]]   #[which.primalchange]
  ord = order(final.model)

  final.model.cluttered = final.model[ord]
  final.model.signs = final.model.signs[ord]
  final.model.early = final.model[final.model <= n]
  final.model.late  = final.model[final.model >  n]
  nbhood=2
  final.model.decluttered = c(declutter(final.model.early,nbhood), final.model.late)
  final.models = list(final.model.cluttered,final.model.decluttered)
#  final.models =  lapply(final.models,function(mymodel) mymodel[mymodel<=n])

  pvals = pvals.decluttered = array(NA,dim=c(n,maxsteps))
  results = results.decluttered = list()
  for(kk in 1){#1:2
    cat("\n type",kk,fill=T)
    final.model = final.models[[kk]]
    if(stop.time > 0){
      ntests = length(final.model)
      for(ii in 1:ntests){
        test.knot      = final.model[ii]        
        adj.knot       = final.model
        test.knot.sign = final.model.signs[which(final.model.cluttered == test.knot)]
        if(test.knot > n) next()
        print(test.knot)  

        
        # Check if there has been a change in primal
        bigmodel = unique(c(adj.knot,test.knot)) #states[[ii]]
        smallmodel = c(adj.knot)
        smallmodel = smallmodel[smallmodel!=test.knot]
        tD = cbind(if(all(is.na(bigmodel))) t(D) else t(D)[,-(bigmodel)])
        tD.prev = cbind(if(all(is.na(smallmodel))) t(D) else t(D)[,-(smallmodel)])
        if( rankMatrix(tD) - rankMatrix(tD.prev) == 0 ) {print("no change in primal, moving on"); next()}

        # make the differences based on the _reduced_ D matrix
        D.reduced = makeDmat(n,order=0)
        adj.knot = adj.knot[adj.knot<=n]
        print(adj.knot)

        # Calculate contrast vector
        v =     make.v.tf.fp(test.knot = test.knot,
                            adj.knot  = adj.knot,
                            test.knot.sign = test.knot.sign,#f1$adj.knot,
                            D = D.reduced)#D=D
        
        coord = final.model[ii]
        pval  = pval.fl1d(y0,G,dik=v,sigma,u=u)
        print(pval)

        if(kk==1){
          results[[ii]] = list(pval=pval,coord=coord,v=v)
        } else {
          results.decluttered[[ii]] = list(pval=pval,coord=coord,v=v)
        }
      }
    }    
  }
  


  save(y0, f0, maxsteps, n, modelinfo, bic, D, consec, sigma, G, u, stop.time, results, results.decluttered,final.model.cluttered,final.model.decluttered,
       file = file.path(outputdir,"cgh-nonsparse-results.Rdata"))



