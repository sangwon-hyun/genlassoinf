  # working directory should be [genlassoinf]/code
  source("settings.R")
  source('funs.R')
  source('examples/testfuns.R')
  source('dualPathSvd2.R')
  lapply(c("genlasso","pryr"), require, character.only = TRUE)
  #library(polypath)

  # 2004 Hot spot paper
    library(cghFLasso)
    data(CGH)
    summary(CGH)


  # 1st half: estimate the model
    y0 = CGH$GBM.y
    even.ind = 2*(1:(length(y0)/2))
    y1 = y0[even.ind]
    y2 = y0[even.ind-1]
    n1 = length(y1)
    n = n1*2
    maxsteps= 50
    D = rbind(makeDmat(n1,order=0))
    tf.obj = trendfilter(y1,ord=0)
    tf.obj.cv = cv.trendfilter(tf.obj,k=5)
    tf.obj$pathobjs
    tf.obj$lambda # in the dual path
    stop.time = which(tf.obj$lambda == tf.obj.cv$lambda.1se) - 1
    tf.model = tf.obj$beta[,stop.time+1]
    tf.cps = tf.obj$pathobj$B[1:stop.time]
    tf.signs = tf.obj$pathobj$s[1:stop.time]
    
    plot(tf.obj.cv)
    plot(tf.obj,lambda=tf.obj.cv$lambda.1se,main="One standard error rule")
    objects(tf.obj)
    lines(tf.model, col='blue')
    abline(v=tf.obj$pathobj$B[1:stop.time],col='grey')
    points(y2, pch=17)


  # 2nd half: conduct z-tests on adjacent points
    ## Conduct z tests (with their respective noises)
  two.sample.test = function(data, testloc, testdir, model, sigma=NA, type=c("z","t")){
    type = match.arg(type)
    n = length(data)
    model.augmented = sort(c(1,model,n))
    test.ind = which(model.augmented == testloc)
    data1 = data[(model.augmented[test.ind-1]+1):(model.augmented[test.ind])]
    data2 = data[(model.augmented[test.ind]+1):(model.augmented[test.ind+1])]
    n1 = length(data1)
    n2 = length(data2)
    
    # Manually form the z/ttests    
    if(type=='z'){
      z.se = sigma* (1/sqrt(n1)+1/sqrt(n2))
      z.stat = (testdir * (mean(data2) - mean(data1))) / z.se 
      pval = 1-pnorm(z.stat,mean=0,sd=1)
      return(pval)
    } else if( type=='t'){
      stop("Bug not resolved yet!")
      # Equal variance two-sample student t-test, unequal sample size
      t.se = sqrt(var(data2)/n2 + var(data1)/n1)
      t.stat = (testdir * (mean(data2) - mean(data1))) / t.se
      pval = 1-pt(t.stat,n1+n2-2)
      # Sanity check for t-test
      if(testdir==1) alt.dir="greater"
      if(testdir==-1) alt.dir="less"
      
      ttest = t.test(data2, data1, alternative=alt.dir, var.equal=T)
      if(!(all.equal(pval, ttest$p.value)==T )){
       print(pval)
       print(ttest)
      }
      t.test(data2, data1, alternative="greater", var.equal=T)$p.value
      return(pval)
    } else { stop("type of test not written yet")}
  }

  #pvals = unlist( Map( function(a,b){ two.sample.test(y1,a,b,tf.cps,sigma,type='t') }, tf.cps, tf.signs))
  pvals = unlist( Map( function(a,b){ two.sample.test(y2,a,b,tf.cps,sigma,type='z') }, tf.cps, tf.signs))
  names(pvals) = tf.cps*2
  round(pvals,3)
  round(pvals[order(tf.cps)],3)
  xtable(rbind(round(pvals[order(tf.cps)],2)))
  pvals[order(tf.cps)]<0.05/length(pvals)

# Visualize the sample splitted detections
  pch.dat=16
  pcol.dat="lightgrey"
  lwd.est = 2
  plot(y0,pch=pch.dat,col=pcol.dat)
  lines(rep(tf.model,each=2),col='blue',lwd=2)
  lines(tf.model, col='blue')
  abline(v=tf.obj$pathobj$B[1:stop.time],col='grey')
  points(y2, pch=17)

## Fiture work: sample splitting doesn't have to be done with 50/50 cut!


## save data

#  pdf(file=file.path(outputdir,"cgh-samplesplitting.Rdata"))
#    CGH.FL.obj1<-cghFLasso(CGH$GBM.y)
##    plot(CGH.FL.obj1, index=1, type="Lines")
#    fit.2004 = CGH.FL.obj1[[1]]
#    lcol.2004 = "red"
#    lwd.2004 = 2
#    dat.2004 = CGH.FL.obj1[[2]]
#    pcol.2004 = "grey50"
#    pch.2004 = 16
#    plot(dat.2004, col = pcol.2004, pch = pch.2004,axes=F)
#    axis(1);axis(2)
#    lines(fit.2004, col = lcol.2004, lwd = lwd.2004)
#    cp.2004 = which(sapply(1:(length(fit.2004)-1), function(ii) fit.2004[ii] != fit.2004[ii+1] ))
#    abline(v=cp.2004, col = 'lightgrey',lty = 3)
#  graphics.off()

## CP package


## Fit fused lasso paths
#  y0 = CGH$GBM.y
#  n = length(y0)
#  maxsteps= 100
#  D = rbind(getDmat(n,order=0),diag(rep(1,n)))

#  n=100
#  beta0 = c(rep(0,n/2),rep(3,n/2))
#  y0 = beta0 + rnorm(n,0,2) 
#  plot(y0)
# 


### Sample Splitting (for p-values)
