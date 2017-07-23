## Make sure you're working from [dropboxfolder]/code
source("settings.R")
source('helpers.R') # generic helpers
source('funs.R')
source('testfuns.R')
source('dualPathSvd2.R')
library(genlasso)
library(RColorBrewer)
outputdir = "output"
codedir = "."

## First, visualize what is happening
  sigma=1
  n = maxsteps = 60
  beta0 = alternjump.y(returnbeta=T, lev1=1, lev2=3, sigma=sigma, n=n)
  y0    = alternjump.y(returnbeta=F, lev1=1, lev2=3, sigma=sigma, n=n)
  f0    = dualpathSvd2(y0, D=dual1d_Dmat(length(y0)), maxsteps,approx=T)
  
  consec =2
  bic = getbic(y0,f0,sigma,maxsteps)
  ind = c()

  plot(y0)
  stoptime = which.rise(bic,consec,n)
  title(paste("stopped at", stoptime-1))
  lines(f0$beta[,ind[jj]],col='blue')
  abline(v = f0$pathobj$B[1:(stoptime-1)],lty=2,col='lightgrey')
  lines(beta0,col='red')

  plot(bic, ylim = c(min(bic),max(bic)*1.5))
  abline(v = stoptime, col = 'lightgrey', lty=2)
  text(x = ind, y = max(bic)+rnorm(length(dirs),0,max(bic)/10), labels = dirs, srt = 90)

  
  
## obtain unconditional power for forward BIC rule
  nsim = 100000
  ngrain = 5
  sigmalist = seq(from= 0.1, to = 5, length= ngrain)
  consec = 2
  nlist = 3*c(10,40,90) #  n = 30 #maxsteps = n/2
  sigmamaxlist = 3*(1:3)
  lev1 = 0
  lev2 = 1
  for(kk in 1:3){
    n = nlist[kk]
    maxsteps = 15  #sigmamaxlist[kk]*2 # this is because.. computation time sucks

    stoptimes.bic = stoptimes.bic2 = stoptimes.ebic = stoptimes.sbic = stoptimes.aic = array(NA, c(ngrain,nsim))
    pvals.bic = pvals.bic.naive = 
    pvals.bic2 = pvals.bic2.naive = 
    pvals.ebic = pvals.ebic.naive = 
    pvals.sbic = pvals.sbic.naive = 
    verdicts.bic = verdicts.bic.naive = 
    verdicts.bic2 = verdicts.bic2.naive = 
    verdicts.ebic = verdicts.ebic.naive = 
    verdicts.sbic = verdicts.sbic.naive = 
    pvals.aic = pvals.aic.naive = 
    verdicts.aic = verdicts.aic.naive = array(NA,c(ngrain, nsim, n))

    pvals.fixed2 = pvals.fixed3 = verdicts.fixed2 = verdicts.fixed3 = array(NA,c(ngrain, nsim, n))
    
    
    pvals.oracle = verdicts.oracle = array(NA,c(ngrain,nsim,n))
    cat('\n', n, "out of", nlist)
    sigmalist = seq(from= 0.1, to = sigmamaxlist[kk], length= ngrain)
    for(ii in 1:ngrain){
      cat('\n', ii, "out of", ngrain)
      sigma = sigmalist[ii]
      cat('\n')
      for(isim in 1:nsim){
        cat('\r', isim, "out of", nsim)

        # generate data + path
        beta0 = alternjump.y(returnbeta=T, lev1=lev1, lev2=lev2, sigma=sigma, n=n)
        y0    = alternjump.y(returnbeta=F, lev1=lev1, lev2=lev2, sigma=sigma, n=n)
        f0    = dualpathSvd2(y0, dual1d_Dmat(length(y0)), maxsteps,approx=T)

        ### 1. collect information
        # bic
        stoptime.bic = which.rise(getbic(y0,f0,sigma,maxsteps),consec, n) - 1 # internally defining the `stoptime' to be the step of the algorithm where you stop. the stoptime to be plotted is this+1.
        stoptime.bic = pmin(stoptime.bic,n-consec-1)
        if(stoptime.bic > 0){ # when stoptime is zero, no test is conducted.  
        locs.bic = f0$pathobj$B[1:stoptime.bic]
        Gobj    = getGammat(f0,y0,stoptime.bic+consec,"dualpathSvd",'bic',sigma,consec,maxsteps)
        G       = Gobj$Gammat
        u       = Gobj$u
        G.naive = getGammat(f0,y0,stoptime.bic+consec,"dualpathSvd",maxsteps=maxsteps)
        for(test.step in 1:stoptime.bic){
          d = getdvec(f0,y0,test.step,stoptime.bic,type="segment",usage="dualpathSvd",matchstep=F)
          loc = locs.bic[test.step]
          pvals.bic[ii, isim, loc]       <- pval.fl1d(y0, G,       d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE, u)
          pvals.bic.naive[ii, isim, loc] <- pval.fl1d(y0, G.naive, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE)
          verdicts.bic[ii, isim, loc]       <- (pvals.bic[ii, isim, loc] < (0.05/stoptime.bic))
          verdicts.bic.naive[ii, isim, loc] <- (pvals.bic.naive[ii, isim, loc] < (0.05/stoptime.bic))
        }
        }
        
          
        stoptime.bic2 = which.rise(getbic(y0,f0,sigma,maxsteps,fac=2),consec, n) - 1 # internally defining the `stoptime' to be the step of the algorithm where you stop. the stoptime to be plotted is this+1.
        stoptime.bic2 = pmin(stoptime.bic2,n-consec-1)
        if(stoptime.bic2 > 0){ # when stoptime is zero, no test is conducted.  
        locs.bic2 = f0$pathobj$B[1:stoptime.bic2]
        Gobj    = getGammat(f0,y0,stoptime.bic2+consec,"dualpathSvd",'bic',sigma,consec,maxsteps)
        G       = Gobj$Gammat
        u       = Gobj$u
        G.naive = getGammat(f0,y0,stoptime.bic2+consec,"dualpathSvd",maxsteps=maxsteps)
        for(test.step in 1:stoptime.bic2){
          d = getdvec(f0,y0,test.step,stoptime.bic2,type="segment",usage="dualpathSvd",matchstep=F)
          loc = locs.bic2[test.step]
          pvals.bic2[ii, isim, loc]       <- pval.fl1d(y0, G,       d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE, u)
          pvals.bic2.naive[ii, isim, loc] <- pval.fl1d(y0, G.naive, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE)
          verdicts.bic2[ii, isim, loc]       <- (pvals.bic2[ii, isim, loc] < (0.05/stoptime.bic2))
          verdicts.bic2.naive[ii, isim, loc] <- (pvals.bic2.naive[ii, isim, loc] < (0.05/stoptime.bic2))
        }
        }
        
        
        # sbic: strengthened bic
        stoptime.sbic = which.rise(getbic(y0,f0,sigma,maxsteps,strength=1.01),consec, n) - 1 # internally defining the `stoptime' to be the step of the algorithm where you stop. the stoptime to be plotted is this+1.
        stoptime.sbic = pmin(stoptime.sbic,n-consec-1)
        if(stoptime.sbic > 0){ # when stoptime is zero, no test is conducted.  
        locs.sbic = f0$pathobj$B[1:stoptime.sbic]
        Gobj    = getGammat(f0,y0,stoptime.sbic+consec,"dualpathSvd",'bic',sigma,consec,maxsteps)
        G       = Gobj$Gammat
        u       = Gobj$u
        G.naive = getGammat(f0,y0,stoptime.sbic+consec,"dualpathSvd",maxsteps=maxsteps)
        for(test.step in 1:stoptime.sbic){
          d = getdvec(f0,y0,test.step,stoptime.sbic,type="segment",usage="dualpathSvd",matchstep=F)
          loc = locs.sbic[test.step]
          pvals.sbic[ii, isim, loc]       <- pval.fl1d(y0, G,       d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE, u)
          pvals.sbic.naive[ii, isim, loc] <- pval.fl1d(y0, G.naive, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE)
          verdicts.sbic[ii, isim, loc]       <- (pvals.sbic[ii, isim, loc] < (0.05/stoptime.sbic))
          verdicts.sbic.naive[ii, isim, loc] <- (pvals.sbic.naive[ii, isim, loc] < (0.05/stoptime.sbic))
        }
        }
        
         
        # ebic: extended bic
        stoptime.ebic = which.rise(getebic(y0,f0,sigma,maxsteps, fac=.5),consec, n) - 1 # internally defining the `stoptime' to be the step of the algorithm where you stop. the stoptime to be plotted is this+1.
        stoptime.ebic = pmin(stoptime.ebic,n-consec-1)
        if(stoptime.ebic > 0){ # when stoptime is zero, no test is conducted.  
        locs.ebic = f0$pathobj$B[1:stoptime.ebic]
        Gobj    = getGammat(f0,y0,stoptime.ebic+consec,"dualpathSvd",'bic',sigma,consec,maxsteps)
        G       = Gobj$Gammat
        u       = Gobj$u
        G.naive = getGammat(f0,y0,stoptime.ebic+consec,"dualpathSvd",maxsteps=maxsteps)
        for(test.step in 1:stoptime.ebic){
          d = getdvec(f0,y0,test.step,stoptime.ebic,type="segment",usage="dualpathSvd",matchstep=F)
          loc = locs.ebic[test.step]
          pvals.ebic[ii, isim, loc]       <- pval.fl1d(y0, G,       d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE, u)
          pvals.ebic.naive[ii, isim, loc] <- pval.fl1d(y0, G.naive, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE)
          verdicts.ebic[ii, isim, loc]       <- (pvals.ebic[ii, isim, loc] < (0.05/stoptime.ebic))
          verdicts.ebic.naive[ii, isim, loc] <- (pvals.ebic.naive[ii, isim, loc] < (0.05/stoptime.ebic))
        }
        }
        
        # aic
        stoptime.aic = which.rise(getaic(y0,f0,sigma,maxsteps),consec, n)-1
        stoptime.aic = pmin(stoptime.aic, n-consec-1)
        if(stoptime.aic > 0){ # when stoptime is zero, no test is conducted.  
        locs.aic = f0$pathobj$B[1:stoptime.aic]
        Gobj    = getGammat(f0,y0,stoptime.aic+consec,"dualpathSvd",'aic',sigma,consec,maxsteps)
        G       = Gobj$Gammat
        u       = Gobj$u
        G.naive = getGammat(f0,y0,stoptime.aic+consec,"dualpathSvd",maxsteps=maxsteps)
        for(test.step in 1:stoptime.aic){
          d = getdvec(f0,y0,test.step,stoptime.aic,type="segment",usage="dualpathSvd",matchstep=F)
          loc = locs.aic[test.step]
          pvals.aic[ii, isim, loc]       <- pval.fl1d(y0, G,       d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE,u)
          pvals.aic.naive[ii, isim, loc] <- pval.fl1d(y0, G.naive, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE)
          verdicts.aic[ii, isim, loc]       <- (pvals.aic[ii, isim, loc]< (0.05/stoptime.aic))
          verdicts.aic.naive[ii, isim, loc] <- (pvals.aic.naive[ii, isim, loc]< (0.05/stoptime.aic))
        }
        }
        
        # fixed stop times
        fixedstoptime = 2
        G.truth = getGammat(f0,y0,fixedstoptime,"dualpathSvd",maxsteps=maxsteps)
        for(test.step in 1:fixedstoptime){
           d = getdvec(f0,y0,test.step,fixedstoptime,type="segment",usage="dualpathSvd",matchstep=F)
           loc = f0$pathobj$B[test.step]
           pvals.fixed2[ii, isim, loc]       <- pval.fl1d(y0, G.truth, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE)
           verdicts.fixed2[ii, isim, loc]       <- (pvals.fixed2[ii, isim, loc] < (0.05/fixedstoptime))
        }
        fixedstoptime = 3
        G.truth = getGammat(f0,y0,fixedstoptime,"dualpathSvd",maxsteps=maxsteps)
        for(test.step in 1:fixedstoptime){
           d = getdvec(f0,y0,test.step,fixedstoptime,type="segment",usage="dualpathSvd",matchstep=F)
           loc = f0$pathobj$B[test.step]
           pvals.fixed3[ii, isim, loc]       <- pval.fl1d(y0, G.truth,  d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE)
           verdicts.fixed3[ii, isim, loc]       <- (pvals.fixed3[ii, isim, loc] < (0.05/fixedstoptime))
        }


        # average stoptimes
        stoptimes.bic[ii,isim] = stoptime.bic
        stoptimes.bic2[ii,isim] = stoptime.bic2
        stoptimes.ebic[ii,isim] = stoptime.ebic
        stoptimes.sbic[ii,isim] = stoptime.sbic
        stoptimes.aic[ii,isim] = stoptime.aic
        
        # oracle
        brks = c(0,n/3,2*n/3,n)
        for(brk.i in 1:(length(brks)-2)){
          loc = brks[brk.i+1]
          dif = abs(mean(y0[(brks[brk.i+1]+1):(brks[brk.i+2])]) - mean(y0[(brks[brk.i]+1):(brks[brk.i+1])]))
          lvl = 0.05/(length(brks)-2)
          n1 = n2 = n/(length(brks)-1)
          z_crit = qnorm(1-lvl)*sigma*sqrt(1/n1 + 1/n2)
          verdicts.oracle[ii,isim,loc] = dif > z_crit
          pvals.oracle[ii,isim,loc]    = 1 - pnorm(dif, mean=0, sd = sigma*sqrt(1/n1^2 + 1/n2^2))
        }
      }
 
        # save
       obj.list1 = c("stoptimes.bic", "stoptimes.bic2", "stoptimes.ebic","stoptimes.sbic", "stoptimes.aic", 
                     "pvals.bic", "pvals.bic.naive", "verdicts.bic", "verdicts.bic.naive",
                     "pvals.bic2", "pvals.bic2.naive", "verdicts.bic2", "verdicts.bic2.naive", 
                     "pvals.ebic", "pvals.ebic.naive", "verdicts.ebic", "verdicts.ebic.naive",  
                     "pvals.sbic", "pvals.sbic.naive", "verdicts.sbic", "verdicts.sbic.naive",  
                     "pvals.aic", "pvals.aic.naive", "verdicts.aic", "verdicts.aic.naive", 
                     "pvals.fixed2", "pvals.fixed3","verdicts.fixed2", "verdicts.fixed3",
                     "pvals.oracle","verdicts.oracle",
                     "sigma","nsim","ngrain","lev1","lev2","n", "sigmalist")
      save(list=obj.list1, file=file.path(outputdir, paste0("bic-alternjump-segmentsize-allbics-temp", n/3, ".Rdata")))
    }
  }
  

# outputdir = "/media/shyun/Bridge/Dropbox/CMU/courses(CURRENT)/genlassoinf/code/maxoutput"
# calculate condit power at each location
  # For each location,
  nlist = 3*c(10,40,90) #  n = 30 #maxsteps = n/2
  for(n in nlist[1]){
    cat('\n', n, "out of", nlist)
    load(file=file.path(outputdir, paste0("bic-alternjump-segmentsize-allbics", n/3, ".Rdata")))

    powers.bic = powers.bic.naive = 
    powers.bic2 = powers.bic2.naive = 
    powers.ebic = powers.ebic.naive = 
    powers.sbic = powers.sbic.naive = 
    powers.aic = powers.aic.naive = array(NA,c(ngrain,n))
    powers.fixed2 = powers.fixed3 = array(NA,c(ngrain,n))
    powers.oracle = powers.proxbic = powers.proxbic2 = powers.proxebic = powers.proxsbic = powers.proxaic = array(NA, ngrain)
    getpow = function(verdicts,ii,loc){  return(sum(verdicts[ii,,loc],na.rm=T)/pmax(1,sum(!is.na(verdicts[ii,,loc]))))   }
    # For each location,
    for(ii in 1:ngrain){
      # exact powers
      for(loc in 1:n){
        powers.bic[ii,loc]        = getpow(verdicts.bic,ii,loc)
        powers.bic.naive[ii,loc]  = getpow(verdicts.bic.naive,ii,loc)
        powers.bic2[ii,loc]       = getpow(verdicts.bic2,ii,loc)
        powers.bic2.naive[ii,loc] = getpow(verdicts.bic2.naive,ii,loc)
        powers.ebic[ii,loc]       = getpow(verdicts.ebic,ii,loc)
        powers.ebic.naive[ii,loc] = getpow(verdicts.ebic.naive,ii,loc)
        powers.sbic[ii,loc]       = getpow(verdicts.sbic,ii,loc)
        powers.sbic.naive[ii,loc] = getpow(verdicts.sbic.naive,ii,loc)
        powers.aic[ii,loc]        = getpow(verdicts.aic,ii,loc)
        powers.aic.naive[ii,loc]  = getpow(verdicts.aic.naive,ii,loc)
        powers.fixed2[ii,loc]     = getpow(verdicts.fixed2, ii,loc)
        powers.fixed3[ii,loc]     = getpow(verdicts.fixed3, ii,loc)  
      }
      
      # approximate powers at a particular true break coordinate
      loc = n/3
      proxwidth = .15*n#log(n) #2
      proxlocs = (loc-proxwidth):(loc+proxwidth)
      proxbic.verdict = proxbic2.verdict = proxebic.verdict = proxsbic.verdict = proxaic.verdict = c()

      for(isim in 1:nsim){
      # bic
        verdicts = verdicts.bic[ii,isim,proxlocs]
        proxbic.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
       # bic2
        verdicts = verdicts.bic2[ii,isim,proxlocs]
        proxbic2.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
       # ebic
        verdicts = verdicts.ebic[ii,isim,proxlocs]
        proxebic.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
       # sbic
        verdicts = verdicts.sbic[ii,isim,proxlocs]
        proxsbic.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
       # aic
        verdicts = verdicts.aic[ii,isim,proxlocs]
        proxaic.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
      }
      powers.proxbic[ii] = sum(proxbic.verdict, na.rm=T)/pmax(sum(!is.na(proxbic.verdict)),1)
      powers.proxbic2[ii] = sum(proxbic2.verdict, na.rm=T)/pmax(sum(!is.na(proxbic2.verdict)),1)
      powers.proxebic[ii] = sum(proxebic.verdict, na.rm=T)/pmax(sum(!is.na(proxebic.verdict)),1)
      powers.proxsbic[ii] = sum(proxsbic.verdict, na.rm=T)/pmax(sum(!is.na(proxsbic.verdict)),1)
      powers.proxaic[ii] = sum(proxaic.verdict, na.rm=T)/pmax(sum(!is.na(proxbic.verdict)),1)

      # oracle power at true break coordinate
      powers.oracle[ii] = sum(verdicts.oracle[ii,,loc],na.rm=T)/nsim
    }
    obj.list2 = c("powers.bic", "powers.bic.naive", 
                  "powers.bic2", "powers.bic2.naive",
                  "powers.ebic", "powers.ebic.naive",
                  "powers.sbic", "powers.sbic.naive",
                  "powers.aic", "powers.aic.naive",
                  "powers.fixed2", "powers.fixed3", "powers.oracle",
                  "powers.proxbic", "powers.proxbic2","powers.proxaic")  
    save(list = c(obj.list1, obj.list2), file=file.path(outputdir, paste0("bic-alternjump-segmentsize-allbics", n/3, ".Rdata")))
  }
  
  
  
  # change stoptimes to conditional values (only the values that were tested)
  nlist = 3*c(10,40,90) #n = 10  
  for(n in nlist[1]){
    loc = n/3
    load(file=file.path(outputdir, paste0("bic-alternjump-segmentsize-allbics", n/3, ".Rdata")))
    stoptimes.bic.cond = stoptimes.bic
    for(jj in 1:ngrain){ stoptimes.bic.cond[jj,is.na(verdicts.bic[jj,,loc])] = NA }
    stoptimes.bic2.cond = stoptimes.bic2
    for(jj in 1:ngrain){ stoptimes.bic2.cond[jj,is.na(verdicts.bic2[jj,,loc])] = NA }
    stoptimes.ebic.cond = stoptimes.ebic
    for(jj in 1:ngrain){ stoptimes.ebic.cond[jj,is.na(verdicts.ebic[jj,,loc])] = NA }
    stoptimes.sbic.cond = stoptimes.sbic
    for(jj in 1:ngrain){ stoptimes.sbic.cond[jj,is.na(verdicts.sbic[jj,,loc])] = NA }
    stoptimes.aic.cond = stoptimes.aic
    for(jj in 1:ngrain){ stoptimes.aic.cond[jj,is.na(verdicts.aic[jj,,loc])] = NA }
   
    obj.list3 = c("stoptimes.bic.cond", "stoptimes.bic2.cond", "stoptimes.ebic.cond", "stoptimes.sbic.cond", "stoptimes.aic.cond")
    save(list=c(obj.list1, obj.list2, obj.list3), file=file.path(outputdir, paste0("bic-alternjump-segmentsize-allbics", n/3, ".Rdata")))
  }
  


## make plot of powers
  nlist = 3*c(10,40,90) #n = 10  
  for(n in nlist[1]){
   cat('\n', n, "out of", nlist)
    load(file=file.path(outputdir, paste0("bic-alternjump-segmentsize-allbics", n/3, ".Rdata")))
    pdf(file.path(outputdir,paste0("bic-alternjump-segmentsize-allbics", n/3, ".pdf")), width=15, height=10)
    
  # plot powers
    loc = n/3
    xlim = c(0,max(sigmalist))
    plot(powers.bic[,loc] ~ sigmalist, type = 'l', lwd=2, ylim = c(0,1), xlim = xlim, axes = F, xlab="noise(sd)", ylab = "condit. powers")
    title("Conditional Power at correct location (n/3), for alternating-jump")
    axis(1); axis(2)

    # plot bic
    lines(powers.bic.naive[,loc] ~ sigmalist, type = 'l', lwd=1)
    lines(powers.proxbic ~ sigmalist, type = 'l', lty=2)
       
    # more lines
    lines(powers.bic2[,loc] ~ sigmalist, type = 'b', lwd=1, pch = "2")
    lines(powers.ebic[,loc] ~ sigmalist, type = 'b', lwd=1, pch = "e")
    lines(powers.sbic[,loc] ~ sigmalist, type = 'b', lwd=1, pch = "s")
    #lines(powers.bic.naive[,loc] ~ sigmalist, type = 'l', lwd=1)
    lines(powers.proxbic ~ sigmalist, type = 'l', lty=2)
    lines(powers.proxbic2 ~ sigmalist, type = 'b', lty=2,pch = "2")
    lines(powers.proxebic ~ sigmalist, type = 'b', lty=2,pch = "e")
    
    # plot AIC
    lines(powers.aic[,loc] ~ sigmalist, type = 'l', col = 'red', lwd=2)
    #lines(powers.aic.naive[,loc] ~ sigmalist, type = 'l', col='red',lwd=1)
    lines(powers.proxaic ~ sigmalist, type = 'l', col = 'red', lty=2)
    
    # plot fixed
    lines(powers.fixed3[,loc] ~ sigmalist, type = 'l', col = 'green', lwd=2)
    lines(powers.fixed2[,loc] ~ sigmalist, type = 'l', col='darkgreen',lwd=2)

  # plot oracle
    lines(powers.oracle~sigmalist, col = 'blue', lwd=2)

      
    # plot bic stoptimes
     addstoptimes = function(stoptimes, sigmalist, col, pch, adjust){
        par(new=TRUE)
        mn = apply(stoptimes,1,mean,na.rm=T) + 1 # adding one because of the way we define stoptime.
        sdev = apply(stoptimes,1,sd,na.rm=T)
        sigmalist_adjusted = sigmalist+rep(adjust,length(sigmalist))
        plot(mn ~ sigmalist_adjusted, pch=pch, ylim = c(0,5), xlim = xlim, axes=F, col = col, xlab = "", ylab = "")
        arrows(sigmalist_adjusted, mn-sdev, sigmalist_adjusted, mn+sdev, length=0.05, angle=90, code=3, col = col)

      }

      addstoptimes(stoptimes.bic,sigmalist, 'lightgrey', 19, 0)
      addstoptimes(stoptimes.bic2,sigmalist, 'lightgrey', "2", 0.05)
      addstoptimes(stoptimes.ebic,sigmalist, 'lightgrey', "e", 0.1)
      addstoptimes(stoptimes.aic,sigmalist, 'darkgrey', 19, -0.05)
      mtext("stoptimes", side = 4, padj = 2)
      axis(4)

    abline(h=c(3,4),lty = c(1,2), col = 'yellow')

  # make legend
    legend("topright", lty = c(rep(c(1,1,2),2),1), lwd = c(rep(c(1,2,1),2),1), legend = c("BIC-naive", "BIC", "BICprox", "AIC-naive", "AIC", "AICprox", "oracle","fixed3", "fixed2"), col = c(rep("black",3), rep("red",3), "blue", "green", "darkgreen"))
    legend("bottomright", pch = c(19,19), col = c("lightgrey", "darkgrey"), legend = c("avg bic stoptime+-1sdv","aic-stoptime"))

  dev.off()
  }
  
  
  

  
  
# plot distribution of stoptimes
  source("settings.R")
  pdf(file.path(outputdir, "stoptimeplots-alternjump-allbic", n/3, ".pdf"), width = 10, height=6)
  par(mfrow = c(2,3))
  library(vioplot)
  #nlist = 2*c(10,40,90) #n = 10  
  for(n in nlist[1]){
    source("settings.R")
    load(file=file.path(outputdir, paste0("bic-onejump-segmentsize-allbics", n/2, ".Rdata")))
    source("settings.R")
    ssl.bic = lapply(1:length(sigmalist), function(ii){stoptimes.bic.cond[ii,] + 1 + rnorm(nsim,0,0.01)})
    ssl.bic = lapply(ssl.bic, function(vec){vec[!is.na(vec)]})
    ssl.bic2 = lapply(1:length(sigmalist), function(ii){stoptimes.bic2.cond[ii,] + 1 + rnorm(nsim,0,0.01)})
        ssl.bic2 = lapply(ssl.bic2, function(vec){vec[!is.na(vec)]})
    ssl.ebic = lapply(1:length(sigmalist), function(ii){stoptimes.ebic.cond[ii,] + 1 + rnorm(nsim,0,0.01)})
        ssl.ebic = lapply(ssl.ebic, function(vec){vec[!is.na(vec)]})
    ssl.sbic = lapply(1:length(sigmalist), function(ii){stoptimes.sbic.cond[ii,] + 1 + rnorm(nsim,0,0.01)})
        ssl.sbic = lapply(ssl.sbic, function(vec){vec[!is.na(vec)]})
    ssl.aic = lapply(1:length(sigmalist), function(ii){stoptimes.aic.cond[ii,] + 1 + rnorm(nsim,0,0.01)})
        ssl.aic = lapply(ssl.aic, function(vec){vec[!is.na(vec)]})

    # bic
    plotic = function(sigmalist, ylim, sslobj, title){
      plot(NA, xlim = c(0,length(sigmalist)+1),ylim=ylim,axes=F, ylab = "stoptimes (1 is null model)", xlab = expression(sigma))
      #title(paste0("bic, n=", n))
      abline(h = 1:length(sigmalist), col = 'lightgrey', lty = 2)
#      sslobj = ssl.ebic
      vioplot2 = function(...){vioplot(...,add=T,na.rm=T)}
      do.call(vioplot2,sslobj)
      axis(side=2); axis(side=1, at = 1:length(sigmalist), labels = signif(sigmalist,2)); title(title)
    }
    plotic(sigmalist, c(0,10), ssl.bic, "bic")
    plotic(sigmalist, c(0,10), ssl.bic2, "bic2")
    plotic(sigmalist, c(0,10), ssl.ebic, "ebic")
    plotic(sigmalist, c(0,10), ssl.sbic, "sbic")
    plotic(sigmalist, c(0,10), ssl.aic, "aic")
  }
  dev.off()

  
# overlay JUST the BIC curves (and the two oracles)
#  pdf(file.path(outputdir, "bic-overlayed-alternjump-fewersim.pdf"), width = 8, height=8)
  pdf(file.path(outputdir, "bic-overlayed-alternjump.pdf"), width = 8, height=8)
  nlist = 3*c(10,40,90) #n = 10  
  for(zz in 1:2){
    zzold = zz
    n = nlist[zz]
    cat('\n', n, "out of", nlist)
    source("settings.R")
    load(file=file.path(outputdir, paste0("bic-alternjump-segmentsize", n/3, ".Rdata")))
    source("settings.R")
    #pdf(file.path(outputdir,paste0("bic-alternjump-segmentsize", n/3, ".pdf")), width=12, height=7)
    zz = zzold  
    # plot powers
    loc = n/3
    xlim = c(0,max(sigmalist))
    print(zz)
    plot(powers.bic[,loc] ~ sigmalist, type = 'l', lwd=2, ylim = c(0,1), xlim = xlim, axes = F, xlab="noise(sd)", ylab = "condit. powers", lty = zz)
    title("Conditional Power at correct location (n/3), for alternating-jump")
    axis(1, padj = zz)
    axis(2)

    # plot oracle
    lines(powers.oracle~sigmalist, col = 'blue', lwd=2, lty = zz)
    par(new=T)
  }

  # make legend
  legend("topright", lty = c(1:3), legend = nlist/3)
  dev.off()
  

# overlay JUST the BIC curves (and the two oracles) on the SAME AXES
  pdf(file.path(outputdir, "bic-overlayed-alternjump-sameaxis.pdf"), width = 8, height=8)
  plot(NA, type = 'l', lwd=2, ylim = c(0,1), xlim = c(0,6), axes = F, xlab="noise(sd)", ylab = "condit. powers", lty = zz)
  axis(1)
  axis(2)
  nlist = 3*c(10,40,90) #n = 10  
  
  for(zz in 1:2){
    zzold = zz
    n = nlist[zz]
    cat('\n', n, "out of", nlist)
    source("settings.R")
    load(file=file.path(outputdir, paste0("bic-alternjump-segmentsize", n/3, ".Rdata")))
    source("settings.R")
    #pdf(file.path(outputdir,paste0("bic-alternjump-segmentsize", n/3, ".pdf")), width=12, height=7)
    zz = zzold  
    # plot powers
    loc = n/3
    xlim = c(0,max(sigmalist))
    print(zz)
    lines(powers.bic[,loc] ~ sigmalist, type = 'l', lwd=2, lty = zz)
    title("Conditional Power at correct location (n/3), for alternating-jump")
    # plot oracle
    lines(powers.oracle~sigmalist, col = 'blue', lwd=2, lty = zz)
  }

  # make legend
  legend("topright", lty = c(1:3), legend = nlist/3)
  dev.off()
  

### plot data example
  load(file=file.path(outputdir, "bic-alternjump.Rdata"))
  source("settings.R")
  pdf(file.path(outputdir,"bic-alternjump-example2.pdf"), width=12, height=5)
  par(mfrow = c(2,3))
    for(n in c(30,120)){
    sigma = .5
    beta0 = alternjump.y(returnbeta=T, lev1=0, lev2=2, sigma=sigma, n=n)  # this could change
    y0    = alternjump.y(returnbeta=F, lev1=0, lev2=2, sigma=sigma, n=n)
#    f0    = dualpathSvd2(y0, dual1d_Dmat(length(y0)), maxsteps,approx=T)
    plot(y0,xlab="",ylab="");lines(beta0,col='red');title(paste("example data with noise=",sigma))
#    lines(fitted(lm(y0~(getH(n,0))[,c(1,11,21)])),col='green')

    sigma = 1
    beta0 = alternjump.y(returnbeta=T, lev1=0, lev2=2, sigma=sigma, n=n)  # this could change
    y0    = alternjump.y(returnbeta=F, lev1=0, lev2=2, sigma=sigma, n=n)
#    f0    = dualpathSvd2(y0, dual1d_Dmat(length(y0)), maxsteps,approx=T)
    plot(y0,xlab="",ylab="");lines(beta0,col='red');title(paste("example data with noise=",sigma))
#    lines(fitted(lm(y0~(getH(n,0))[,c(1,11,21)])),col='green')
    
    sigma = 2
    beta0 = alternjump.y(returnbeta=T, lev1=0, lev2=2, sigma=sigma, n=n)  # this could change
    y0    = alternjump.y(returnbeta=F, lev1=0, lev2=2, sigma=sigma, n=n)
#    f0    = dualpathSvd2(y0, dual1d_Dmat(length(y0)), maxsteps,approx=T)
    plot(y0,xlab="",ylab="");lines(beta0,col='red');title(paste("example data with noise=",sigma))
#    lines(fitted(lm(y0~(getH(n,0))[,c(1,11,21)])),col='green')
    }
dev.off()




