# working directory should be [genlassoinf]/code
workingdir = '~/../../media/shyun/Bridge/Dropbox/CMU/courses(CURRENT)/genlassoinf/code'
setwd(workingdir)
source("settings.R")
source('helpers.R') # generic helpers
source('funs.R')
source('examples/testfuns.R')
#
source('dualPathSvd2.R')
library(genlasso)
library(polypath)

# call data
  library(cghFLasso)
  data(CGH)
  summary(CGH)
  plot(CGH$GBM.y)
  

# Fit fused lasso paths
  y0 = CGH$GBM.y
  n = length(y0)
  maxsteps= 30
  #f0 = dualpathSvd2(y0,dual1d_Dmat(n),maxsteps=maxsteps,approx=T,verbose=T)
  #alpha = 0.1
#  save(y0, f0, maxsteps, n, file=file.path(outputdir,"cgh.Rdata"))
    
  load(file=file.path(outputdir,"cgh.Rdata"))
  tf.order = 0
  consec=2
  
# Use cv to recover the standard deviation, for our procedure.
  tf.obj = trendfilter(y0,ord=0,maxsteps=30)
  kfold = 10
  tf.obj.cv = cv.trendfilter(tf.obj, k=kfold, mode = "lambda", approx=F, verbose=T)
  min.1se.lambda = tf.obj.cv$lambda.1se #   tf.obj.cv$lambda.min
  plot(tf.obj,lambda=tf.obj.cv$lambda.1se,main="Minimal CV error")
  stop.time = which(tf.obj$lambda == min.1se.lambda) - 1
  lines(tf.obj$beta[,stop.time+1],col='blue')
  RSS = sum((y0 - tf.obj$beta[,stop.time+1])^2)
  sigma = sqrt(RSS/(n-kfold))
  
# Get stopping time and Polyhedron
  bic = getbic(y0, f0, sigma,maxsteps=maxsteps,k=0)
  stop.time = which.rise(bic,consec)-1
  



  
# No decluttering
  signs = f0$pathobjs$s[1:stop.time]
  states = get.states(f0$action)
  final.model = states[[stop.time]]
  #final.model = declutter(final.model, 1)
  signs = signs[order(final.model)]
  signs = c(-signs[1], signs)
  signs.cumsum = cumsum(signs)
  final.model = sort(final.model)
  indices = Map(function(a,b){a:b},final.model[1:(length(final.model)-1)]+1, final.model[2:length(final.model)])
  indices[2:(length(indices)+1)] = indices
  indices[[length(indices)+1]] = (final.model[length(final.model)]+1):length(y0)
  indices[[1]] = 1:final.model[1]
  sign0 = do.call(c, sapply(1:length(signs.cumsum), function(ii){rep(signs.cumsum[ii], length(indices[[ii]]))}))





# Using _liberalized_ declutter by 5 and 10, and consecutive-sign decluttering (TODO: make this into a function)
  how.close = 10
  signs = f0$pathobjs$s[1:stop.time]
  states = get.states(f0$action)
  final.model.before.declutter = states[[stop.time]]
  # Just declutter by closeness
  final.model.declutter.vanilla = declutter(final.model.before.declutter, how.close= how.close, sort=F)

  which.declutter.eliminate = (1:stop.time)[! (1:stop.time) %in% declutter(final.model.before.declutter, how.close= how.close, sort=F, indexonly=T)]    
  # Declutter by closeness and sign repetitiveness
  # function to further check for any sign continuity
  clutters.sign.too = function(index, signs,final.model.before.declutter, how.close){
    eliminate = F
    next.to.index = order(abs(final.model.before.declutter[index] - final.model.before.declutter))[c(2,3)]
    next.to.index[3] = next.to.index[2]
    next.to.index[2] = index
    
    in.range = ( next.to.index>=0 & next.to.index<= stop.time)
    signs.next.to.index = signs[next.to.index]
    signs.next.to.index[!in.range] = NA
    next.to.index[!in.range] = NA
    ## 2 and 1
    if(!is.na(signs.next.to.index[1]) &
      signs.next.to.index[1] == signs.next.to.index[2] & # signs are close
      abs(final.model.before.declutter[next.to.index[1]] - final.model.before.declutter[next.to.index[2]]) <= how.close  ){ # locations are close
      eliminate=T
     }
    ## 2 and 3
    if(!is.na(signs.next.to.index[3]) &
       signs.next.to.index[2] == signs.next.to.index[3] & # signs are close
       abs(final.model.before.declutter[next.to.index[2]] - final.model.before.declutter[next.to.index[3]]) <= how.close  ){ # locations are close
      eliminate=T
    }
    return(eliminate)#return(c(which.declutter.eliminate, further.eliminate))
  }
  to.eliminate = which.declutter.eliminate[sapply(which.declutter.eliminate, function(index) clutters.sign.too(index,signs,final.model.before.declutter,how.close))]
#   clutters.sign.too(10,signs,final.model.before.declutter,how.close)
  final.model.declutter.loc.and.sign = final.model.before.declutter[-to.eliminate]

  sign0.decluttered = list()
  for(ii in 1:2){
    signs = f0$pathobjs$s[1:stop.time]
    final.model = list(final.model.declutter.vanilla, final.model.declutter.loc.and.sign)[[ii]]
    signs = signs[which(final.model.before.declutter %in% final.model)]
    signs = signs[order(final.model)]
    signs = c(-signs[1], signs)

    signs.cumsum = cumsum(signs)
    final.model = sort(final.model)
    indices = Map(function(a,b){a:b},final.model[1:(length(final.model)-1)]+1, final.model[2:length(final.model)])
    indices[2:(length(indices)+1)] = indices
    indices[[length(indices)+1]] = (final.model[length(final.model)]+1):length(y0)
    indices[[1]] = 1:final.model[1]
    sign0.decluttered[[ii]] = do.call(c, sapply(1:length(signs.cumsum), function(ii){rep(signs.cumsum[ii], length(indices[[ii]]))}))
  }



##########################
### Make sign-set plot ###
##########################
  pdf(file=file.path(outputdir, "signpic.pdf"), width=10,height=6)
  par(mfrow = c(3,1))
    # plot data and fitted final model
    par(mar=c(5.1,4.1,3.1,2.1)-1.5)
    plot(y0, ylab = "Data and final model", xlab="", main = "Data and final model", pch=16,cex=.5,axes=F) 
    axis(1); axis(2)
    legend("topright", pch= c(NA,1), pt.cex=c(NA,.5),lty=c(1,NA),col = c("blue","black"),legend=c("fitted final fused lasso model","data"))
    lines(f0$beta[,stop.time+1],col = 'blue')
    abline(v=final.model.before.declutter,col='lightgrey',lty=2)
    
    # plot sign pictures
    par(mar=c(3.1,4.1,3.1,2.1)-1.5)
    plot(sign0,type='l', axes=F, xlab = "Coordinates")
    title(main = "Sign set plot", line = -2)
#    axis(1); axis(2)
    abline(v=final.model.before.declutter,col='lightgrey',lty=2)
#    # plot vanilla de-cluttered picture
#    plot(sign0.decluttered[[1]],type='l', main = "After decluttering (by 10)")
#    abline(v=100*(0:10),col='lightgrey',lty=2)
    # plot loc and sign de-cluttered sign picture
    plot(sign0.decluttered[[2]],type='l', axes=F)
    title(main = "After decluttering (by neighborhood of 10) \n with sign filtering", line=-2)
    abline(v=final.model.declutter.loc.and.sign,col='lightgrey',lty=2)
  dev.off()



# Do inference
  plot(log(bic)); abline(v=stop.time+1, lty=2)
  bic.fac=1
  Gobj.new.with.stoptime = getGammat.with.stoprule(obj=f0,y=y0,
                                     condition.step = stop.time + consec,
                                     stoprule = "bic", sigma=sigma, bic.fac=bic.fac, type = "tf",
                                     consec=consec, maxsteps=maxsteps, tf.order=tf.order, usage="dualpathSvd")

  G = G.bic = Gobj.new.with.stoptime$Gammat
  u = u.bic = Gobj.new.with.stoptime$u
  
#  G.naive = getGammat.naive(f0,y0,4)$G
  
        # Check correctness of polyhedron:
  #      print(tail(sort(G%*%y0 - u, decreasing=T)))
  #      print(tail(G%*%y0 - u))
  #      all((G%*%y0 > u))

  # Conduct tests and record pvals + stoptimes
  states = get.states(f0$action)
        
  # test pre / post decluttered states, with /their/ adjusted contrasts
  final.models = list(final.model.before.declutter,final.model.declutter.vanilla,   final.model.declutter.loc.and.sign)
  myresult.bic = myresult.bic.decluttered.vanilla= myresult.bic.decluttered.loc.and.sign = array(NA,dim=c(n,nsim,2))
  H = getH.trendfilter(n,tf.order)
  for(kk in 1:3){
    print(kk)
    final.model = final.models[[kk]]
    if(stop.time > 0){
      for(ii in 1:length(final.model)){
      for(ll in 1:2){
      print(ii)
        test.coord = final.model[ii]+1
        adj.coord = c(1:(tf.order+1), (1+final.model))
        adj.coord = adj.coord[adj.coord!=test.coord]
        X = H[, adj.coord ]
#              PH1 = X %*% solve(t(X) %*% X , t(X) )
        In = diag(length(y0))
        PH = fitted(lm(In ~ X -1))
        orthPH = (diag(n) - PH)
        
        # Get Segment P-value
        v = as.numeric(H[,test.coord] %*% orthPH)
        if( v %*% y0 < 0) v = -v
        coord = final.model[ii]
        pval.segment  = pval.fl1d(y0,G,dik=v,sigma,u=u)
        #(pval.segment.naive = pval.fl1d(y0,G.naive, dik=v, sigma=sigma, u = rep(0,nrow(G.naive))))

        # Get Spike P-value
        v.spike=v
        v.spike[abs(v)<1E-10]=0
        cp = which(abs(v.spike[1:(length(v.spike)-1)]-v.spike[2:length(v.spike)]) > 1E-5)[2]
        v.spike[cp]=1; v.spike[cp+1]=-1
        v.spike[-c(cp,cp+1)]=0
        if(v.spike%*%y0 < 0) v.spike=-v.spike
        as.numeric(H[,test.coord] %*% orthPH)
        if( v %*% y0 < 0) v = -v
        
        pval.spike  = pval.fl1d(y0,G,dik=v.spike,sigma,u=u)

        pvals = c(pval.segment,pval.spike)
        pval = pvals[ll]

        if(kk==1){
          myresult.bic[coord,ii,ll] = pval
        } else if (kk==2){
          myresult.bic.decluttered.vanilla[coord,ii,ll] = pval
        } else {
          myresult.bic.decluttered.loc.and.sign[coord,ii,ll] = pval
        }
      }}
    }
  }
  myresults = list(myresult.bic,myresult.bic.decluttered.vanilla, myresult.bic.decluttered.loc.and.sign)

pdf(file=file.path(outputdir, "cgh-with-decluttering.pdf"),width=10,height=10)
  par(mfrow=c(3,1))
  for(kk in 1:3){
    print(kk)
   
    final.model = final.models[[kk]]
    myresult = myresults[[kk]]
    signs = (if(kk==1) sign0 else if (kk%in%c(2,3)) sign0.decluttered[[kk-1]])
    aa.segment = apply(myresult[,,1], 1, function(row){ if(all(is.na(row))) NA else row[!is.na(row)] })    
    aa.spike = apply(myresult[,,2], 1, function(row){ if(all(is.na(row))) NA else row[!is.na(row)] })    

    plot(NA,xlim=c(0,n), cex=.2,xlab='coordinate',ylab='y', ylim = c(-7,5))#ylim=range(y0)
    abline(v=f0$pathobj$B[1:(stop.time)],col='lightgrey',lty=3)
    points(y0, cex=.2)
    lines(f0$beta[,stop.time+1],col='red',lwd=2)
    legend("topright", lty = 1, col='skyblue', legend = "fl-detected jump",bg='white')
        
    # overlay sign plot
    lines(signs, col='cyan')
    
#    text(x=final.model, y = 3, labels = paste('seg test p-value:\n',pvals.segment[algstep]))
    text(x = final.model, y = signs[final.model+1]+.5-1, labels = round(aa.segment[final.model],3))#rep(-3, length(final.model))+runif(length(final.model),-1,1)
    text(x = final.model, y = signs[final.model+1]-1, labels = round(aa.spike[final.model],3), col = 'darkblue')#rep(-3, length(final.model))+runif(length(final.model),-1,1)
    text(x = sort(final.models[[3]])[c(1:4,8)], y = -6.3, labels = c("A","B","C","D","E"))
  }
dev.off()
  

## Single oracle testfuns
plot(c(y0[251:737],y0[738:990]))
abline(v = 737.5-250)
mysd = sqrt((1/length(738:990) + 1/length(251:737)) * sigma^2)
myZ = (mean(y0[738:990]) - mean(y0[251:737]) )
oracle.pval = 1- pnorm(myZ, mean=0, sd = mysd)
signif(oracle.pval,10)

# 


  
 
#  

#    text(x=f0$pathobj$B[algstep], y = -2, labels = paste('spike test p-value:\n',pvals.spike[algstep]),)
#          
#    # Calculate p-values
#    pvals.spike[algstep]   = signif(pval.fl1d(y=y0, G=G, dik=d.spike, sigma=sigma, approx=TRUE, threshold=TRUE, approxtype="rob"),2)
#    pvals.segment[algstep] = signif(pval.fl1d(y=y0, G=G, dik=d.segment, sigma=sigma, approx=TRUE, threshold=TRUE, approxtype="rob"),2)
#    if(pvals.spike[algstep] < alpha) spike.sigsteps = c(spike.sigsteps, step)
#    if(pvals.segment[algstep] < alpha) seg.sigsteps = c(seg.sigsteps, step)
#    

#  save(file=file.path(outputdir,"cgh.Rdata"), list = ls())



## Plot some spike results
#  load(file=file.path(outputdir,"cgh.Rdata"))
#  pdf(file.path(outputdir,"cgh-spike-afew.pdf"),width=15,height=4)
#  par(mfrow=c(1,4))
#  for(algstep in spike.steplist){#
#    step = algstep+1
#    plot(NA,xlim=c(0,n),ylim=range(y0), cex=.2,xlab='coordinate',ylab='y')
#    abline(v=path$pathobj$B[1:(algstep-1)],col='lightgrey',lty=3)
#    abline(v=path$pathobj$B[algstep],col='skyblue')
##    abline(v=path$pathobj$B[sigsteps],col='darkseagreen1')
#    points(y0, cex=.2)
#    lines(path$beta[,step],col='red',lwd=2)
#    legend("topright", lty = 1, col='skyblue', legend = "fl-detected jump",bg='white')
#    text(x=path$pathobj$B[algstep], y = 3, labels = paste('seg test p-value:\n',pvals.segment[algstep]))
#    text(x=path$pathobj$B[algstep], y = -2, labels = paste('spike test p-value:\n',pvals.spike[algstep]))
#    # title(paste0("spike pvalue: ", pvals.spike[algstep], "\n", "segment pvalue:",pvals.segment[algstep]))
#  }
#  dev.off()

## plot some segment results
#  load(file=file.path(outputdir,"cgh.Rdata"))
#  pdf(file.path(outputdir,"cgh-segment-afew.pdf"),width=15,height=4)
#  par(mfrow=c(1,4))
#  for(algstep in seg.steplist){
#    step = algstep+1
#    plot(NA,xlim=c(0,n),ylim=range(y0), cex=.2,xlab='coordinate',ylab='y')
#    abline(v=path$pathobj$B[1:(algstep-1)],col='lightgrey',lty=3)
#    abline(v=path$pathobj$B[algstep],col='skyblue')
#    points(y0, cex=.2)
#    lines(path$beta[,step],col='red',lwd=2)
#    legend("topright", lty = 1, col='skyblue', legend = "fl-detected jump",bg='white')
#    text(x=path$pathobj$B[algstep], y = 3, labels = paste('seg test p-value:\n',pvals.segment[algstep]))
#    text(x=path$pathobj$B[algstep], y = -2, labels = paste('spike test p-value:\n',pvals.spike[algstep]))
#    # title(paste0("spike pvalue: ", pvals.spike[algstep], "\n", "segment pvalue:",pvals.segment[algstep]))
#  }
#  dev.off()
#  
## plotting the results after about 15 steps, side by side
#  load(file=file.path(outputdir,"cgh.Rdata"))
#  pdf(file.path(outputdir,"cgh-compare.pdf"),width=14,height=7)
#    par(mfrow=c(2,1))
#    plot(NA,xlim=c(0,n),ylim=range(y0), cex=.2,xlab='coordinate',ylab='y')
#    abline(v=path$pathobj$B[spike.sigsteps],col='blue',lty=2)
#    points(y0, cex=.2)
#    lines(path$beta[,algstep],col='red',lwd=2)
#    legend("topright", lty = 1, col='blue', legend = "spike signif",bg='white')

#    plot(NA,xlim=c(0,n),ylim=range(y0), cex=.2,xlab='coordinate',ylab='y')
#    abline(v=path$pathobj$B[seg.sigsteps],col='pink',lty=2)
#    points(y0, cex=.2)
#    lines(path$beta[,step],col='red',lwd=2)
#    legend("topright", lty = 1, col='pink', legend = "segment signif",bg='white')
#  dev.off()
#  
#  
#  
#  
#  
### What happens when we condition on a further step, for segment?
## generate and save output
#  pvals.spike.further=pvals.segment.further=c()
##  seg.steplist = c(1,4,5,10)
##  spike.steplist=c(5,7,9,10)
#  steplist=1:16#seg.steplist#spike.steplist#c(1:16)
#  seg.sigsteps = spike.sigsteps = c()
#  for(algstep in steplist){
#    step = algstep+1
#    sigma = sd(y0-path$beta[,step])
#    d.spike   = getdvec(obj=path, y=y0, k=algstep, usage = "dualpathSvd", type="spike", matchstep=TRUE)
#    d.segment = getdvec(obj=path, y=y0, k=algstep, klater = step, usage = "dualpathSvd", type="segment", matchstep=F)
#    G = path$Gammat[1:path$nk[algstep+1],]
#    pvals.spike[algstep]   = signif(pval.fl1d(y=y0, G=G, dik=d.spike, sigma=sigma, approx=TRUE, threshold=TRUE, approxtype="rob"),2)
#    pvals.segment.further[algstep] = signif(pval.fl1d(y=y0, G=G, dik=d.segment, sigma=sigma, approx=TRUE, threshold=TRUE, approxtype="rob"),2)
#    if(pvals.spike[algstep] < alpha) spike.sigsteps = c(spike.sigsteps, step)
#    if(pvals.segment[algstep] < alpha) seg.sigsteps = c(seg.sigsteps, step)
#  }
#  save(file=file.path(outputdir,"cgh.Rdata"), list = ls())



## comparing the difference (doesn't do well; understandably so.)
#  pdf(file.path(outputdir,"cgh-segmentcompare.pdf"),width=15,height=4)
#  par(mfrow=c(1,4))
#  for(algstep in seg.steplist){
#    step = algstep+1
#    plot(NA,xlim=c(0,n),ylim=range(y0), cex=.2,xlab='coordinate',ylab='y')
#    abline(v=path$pathobj$B[1:(algstep-1)],col='lightgrey',lty=3)
#    abline(v=path$pathobj$B[algstep],col='skyblue')
#    points(y0, cex=.2)
#    lines(path$beta[,step],col='red',lwd=2)
#    legend("topright", lty = 1, col='skyblue', legend = "fl-detected jump",bg='white')
#    text(x=path$pathobj$B[algstep], y = 3, labels = paste('seg test p-value:\n',pvals.segment[algstep]))
#    text(x=path$pathobj$B[algstep], y = 2.5, labels = paste('1-more-condit p-value:\n',pvals.segment.further[algstep]))
#    # title(paste0("spike pvalue: ", pvals.spike[algstep], "\n", "segment pvalue:",pvals.segment[algstep]))
#  }
#  dev.off()


## what happens when we condition on step 16, and do everything? 

## The key question is now what exactly is step 16?
