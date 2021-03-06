## Make sure you're working from [dropboxfolder]/code
source("settings.R")
source('funs.R')
source('testfuns.R')
source('selectinf/selectiveInference/R/funs.inf.R')
source('dualPathSvd2.R')
library("stepR")
library("Hmisc")
library(cghFLasso)
data(CGH)

outputdir = "output"

## Try it once  
y0 = CGH$GBM.y
sy = smuceR(y0,jumpint=T,confband=T)
cb = confband(sy)
plot(y0,col='grey50')
lines(sy)
lines(cb,col='blue')


## obtain things that will help calculate power
  consec = 2
  n=60
  maxsteps = n-2
  D = makeDmat(n,'trend.filtering',ord=0)
  ngrain = 10
  lev2list = seq(from=0, to=4,length=ngrain)
  sigma = 1
  lev1=0 ; lev2 = 1
  n=60 ;   sigma=1
  nsimlist = seq(from=200,to=4000,length=ngrain)
  
  contrast.smuce = function(y,d, alpha=.05){
    v = d/sqrt(sum(d^2))
    sm = smuceR(y,jumpint=T,confband=T,alpha=alpha)
    cb = confband(sm)
    cb.unordered.pointwise = cb[,2:3] * v
    cb.ordered.pointwise = apply(cb.unordered.pointwise,1,function(cc)range(cc) )
    return(apply(cb.ordered.pointwise,1,sum))
  }
  
  pvals.list = cis.list = list()
  for(igrain in 1:ngrain){
    #sigma = sigmalist[igrain]
    lev2 = lev2list[igrain]
    maxsteps = 5+2*ngrain
    nsim = nsimlist[igrain]
    cat('\n', "lev2 is", lev2, '\n')
    pvals = array(NA,dim=c(nsim,n))  
    cis = array(NA,dim=c(nsim,n,2))

    for(isim in 1:nsim){  
      cat(isim," ")
      beta0 = rep(c(lev1,lev2),each=n/2)
      D = makeDmat(n,'trend.filtering',ord=0)
      y0    = beta0 + rnorm(n, 0,sigma)
      f0    = dualpathSvd2(y0,D,maxsteps,approx=T)

      mm = get.modelinfo(f0,y0,sigma,maxsteps,D=D,stoprule='bic')
      bic = mm$ic      
      stoptime.bic = which.rise(bic,consec) - 1 # internally defining the `stoptime' to be the step of the algorithm where you stop. the stoptime to be plotted is this+1.
      stoptime.bic = pmin(stoptime.bic, n-consec-1)
      if(stoptime.bic == 0) next
          
      locs.bic = f0$pathobj$B[1:stoptime.bic]
      Gobj    = getGammat.with.stoprule(obj=f0,y=y0,condition.step=stoptime.bic+consec,
                                        type='tf',stoprule='bic',sigma=sigma,consec=consec,maxsteps=maxsteps, D=D)
      G       = Gobj$Gammat
      u       = Gobj$u

      # Non-decluttered
      for(test.step in 1:stoptime.bic){
        loc = locs.bic[test.step]
        d       = getdvec(f0,y0,test.step,stoptime.bic,type="segment")
        pval    = pval.fl1d(y0, G, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE, u)
        ci.smuce    = contrast.smuce(y0,d)
        pvals[isim,loc] = pval
        cis[isim,loc,1:2] = ci.smuce
      }
    }
    pvals.list[[igrain]] = pvals
    cis.list[[igrain]] = cis
  }  
  
  save(pvals.list,cis.list,consec,ngrain,lev1,lev2,maxsteps,D,ngrain,sigmalist,nsimlist,n,
       file = file.path(outputdir,"smuce.Rdata"))

########################################################
### PLOT + calculate conditional power #################
########################################################

  load(file=file.path(outputdir,"smuce.Rdata"))
  for(loc in 30:31){
    pow.postsel = c()
    pow.smuce = c()

    for(igrain in 1:ngrain){
      pvals = pvals.list[[igrain]]
      cis = cis.list[[igrain]]
      
      stoptimes = apply(pvals,1,function(myrow)sum(!is.na(myrow)))  
      verdicts = unlist(Map(function(pval,stoptime){ (pval < 0.05/stoptime)}, pvals[,loc], stoptimes))
      
      a = (cis[,loc,1:2])
      verdicts.smuce = apply(a,1,function(myrow){return(0<min(myrow)| 0>max(myrow))})
        
      ## Get powers
      pow.postsel[igrain] = sum(verdicts,na.rm=T)/sum(!is.na(verdicts))
      pow.smuce[igrain]   = sum(verdicts.smuce,na.rm=T)/sum(!is.na(verdicts.smuce))
      if(loc!=30) pow.postsel[1] = pow.smuce[1] = NA
    }

    mar = c(4.5,4.5,2.53,0.5)
    w=5; h=5.5
    pdf(file=file.path(outputdir,paste0("condit-power-smuce-loc",loc,".pdf")), width=w,height=h)
    ## plot power
    par(mar=mar)
    pchs=c(16,17)
    cols=c('black','red')
    legends=c('TG','SMUCE')
    main=bquote(atop(Power, (location==.(loc)))) 
    pows = list(pow.postsel, pow.smuce)
    denom.smuce =  sapply(1:ngrain, function(igrain) sum(!is.na(pvals.list[[igrain]][,loc])))
    denom.fl =  sapply(1:ngrain, function(igrain) sum(!is.na(  cis.list[[igrain]][,loc,1])))  
    denoms = list(denom.smuce,denom.fl)
    lwd.pow=1
    for(ii in 1:length(pows)){
      pow = pows[[ii]]
      if(ii==1){
        plot(pow~lev2list,ylim=c(0,1),type='o',pch=pchs[ii],ylab="Power",xlab=bquote(delta),col=cols[1],axes=F,lwd=lwd.pow)
        title(main=main)
        axis(1); axis(2)
        legend("topleft", lty=c(1,1),pch=pchs,col=cols,legend=legends,lwd=rep(lwd.pow,2))
      } else {
        lines(pow.smuce~lev2list,type='o',pch=pchs[ii],col=cols[ii],lwd=lwd.pow)
      }
      denom = denoms[[ii]]
      binom.std.err = unlist(Map(function(p,n) sqrt(p*(1-p)/n), pow, denom))
      middle = pow
      upper = pow+1.96*binom.std.err
      lower = pow-1.96*binom.std.err
      errbar(x = lev2list[length(lev2list):1], 
             y = middle[length(middle):1], 
             yplus = upper[length(upper):1],
             yminus = lower[length(lower):1],
             pch=pchs[ii],col=cols[ii],errbar.col = cols[ii],add=T)
             

    }
   graphics.off()
  }

############################
#### Plot data  ############
############################
  mar = c(4.5,4.5,2.5,0.5)
  w=5; h=5.5
  pdf(file=file.path(outputdir,paste0("smuce-data.pdf")), width=w,height=h)
  sigma=1
  par(mar=mar)
  set.seed(1)
  n=60
  lev2=2
  beta0 = rep(c(lev1,lev2),each=n/2)
  y0 = beta0 + rnorm(n, 0,sigma)
  D = makeDmat(n,'trend.filtering',ord=0)
  f0    = dualpathSvd2(y0,D,maxsteps,approx=T)
  mm = get.modelinfo(f0,y0,sigma,maxsteps,D=D,stoprule='bic')
  bic = mm$ic   
  consec=2   
  stoptime.bic = which.rise(bic,consec) - 1 
  fl.est = f0$beta[,stoptime.bic+1]
  pcol.dat='grey50'
  pch.dat = 16
  lcol.sig = 'red'
  lwd.sig = 3
  col.fl = 'red'
  col.smuce='blue'
  col.smuce.conf.band = 'grey90'
  lwd.fl=3
  lwd.smuce.conf.band=1
  lwd.smuce.est=3
  col.null.contrast = 'green'
  ylab = ""
  xlab = "Location"
  ylim = c(min(y0),max(y0)*1.3)
  xlim = c(0,70)
  d = c(rep(0,n/4),rep(1,n/4),rep(0,n/2))
  main = bquote(atop(Data~example))
  plot(NA,main=main,axes=F,ylab=ylab,xlab=xlab,ylim=ylim,xlim=xlim)

  
#  matlines(cb[,2:3],col=col.smuce,lwd=lwd.smuce.conf.band,lty=1)
  # shade the area instead  
  sm = smuceR(y0,jumpint=T,confband=T)
  cb = confband(sm)
  y.low <- cb[,2]
  y.high <- cb[,3]
  
  polygon(c(1:nrow(cb), rev(1:nrow(cb))), c(y.high, rev(y.low)),
          col = col.smuce.conf.band, border = NA)

  lines(beta0,col=lcol.sig,lwd=lwd.sig)
  lines(sm,lwd=lwd.smuce.est,col=col.smuce)

  points(y0,pch=pch.dat,col=pcol.dat)
  #lines(fl.est, lwd=lwd.fl,col=col.fl)
  #lines(d,col=col.null.contrast)
  #text(x = c(65,65,65), y = c(0,1,2), labels = c(expression(paste(delta,"=0")),expression(paste(delta,"=1")),expression(paste(delta,"=2"))),cex= c(1.2,1.2,1.2))
  text(x = c(65,65),  y = c(0,lev2), labels = c(expression(paste(delta,"=",0)),expression(paste(delta,"=", 2))),cex= c(1.2,1.2))
  arrows(x0=60-1, y0=0+0.05, x1 = 60-1, y1 = lev2-0.05, code = 3, length = 0.1, angle = 30, lty=2)
  
  legend(x=0,y=5.1, col=c(lcol.sig,col.smuce,col.smuce.conf.band),lwd=c(lwd.sig,lwd.smuce.est,10), legend = c('Signal','SMUCE estimate', 'SMUCE confidence band (95%)'))
  axis(1)
  axis(2)
graphics.off()


#################################################
## Obtain Type 1 error, by n, for SMUCE vs TF ###
#################################################
  cis.null.list = pvals.null.list = list()
  nsim = 500
  sigma = 1
  nlist = c(1:5)*60
  ngrain=length(nlist)
  lev1=0
  lev2=1
  for(igrain in 1:ngrain){
    n = nlist[igrain]
    cat('\n', "n is", n, '\n')
    cis = array(NA,dim=c(nsim,2))
    pvals = rep(NA,nsim)
    D = makeDmat(n,'trend.filtering',ord=0)

    for(isim in 1:nsim){  
      cat("\r", isim, " ")
      beta0 = rep(c(lev1,lev2),each=n/2)
      y0    = beta0 + rnorm(n, 0,sigma)

      ## true contrast
      d = c(rep(-1,n/2), rep(+1,n/2))
      d = d/sqrt(sum(d^2))

      ## false contrast
      d = c(rep(0,n/4),rep(+1,n/4),rep(0,n/2))
      d = d/sqrt(sum(d^2))

      
      ## SMUCE
      ci.smuce    = contrast.smuce(y0, d, .05)  
      cis[isim,1:2] = ci.smuce
      
      ## TG
      maxsteps=3+1
      f0  = dualpathSvd2(y0,D,maxsteps,approx=T)
      stoptime =3
      if(stoptime == 0) next
          
      Gobj    = getGammat.naive(obj=f0,y=y0,condition.step=stoptime)
      G       = Gobj$Gammat
      u       = rep(0,nrow(G))

      pval    = pval.fl1d(y0, G, d, sigma, approx=TRUE, approxtype = "rob", threshold=TRUE, u)
      pvals[isim] = pval
    }
      
    cat(fill=T)
    cis.null.list[[igrain]] = cis
    pvals.null.list[[igrain]] = pvals

    save(pvals.null.list,cis.null.list,ngrain,lev1,lev2,nlist,D,ngrain,sigma,
         lev1,lev2,nsim, file = file.path(outputdir,"smuce-null.Rdata"))
  }


  # Calculate type 1 error  
  load(file = file.path(outputdir,"smuce-null.Rdata"))
  type.1.errors.smuce = type.1.errors.TG = c()
  for(igrain in 1:ngrain){
      type.1.errors.TG[igrain]    = 1 - sum(apply(cis.null.list[[igrain]],1,
                                                  function(myrow) return( min(myrow) <0 & 0 < max(myrow) )))/nsim
    type.1.errors.smuce[igrain] = sum(sapply(pvals.null.list[[igrain]],function(myvec) return( sum(myvec < 0.05) )))/nsim#sum(apply(cis.null.list[[igrain]],1,function(myrow) return( min(myrow) <0 & 0 < max(myrow) )))/nsim
  }
  errs.list = list(type.1.errors.TG, type.1.errors.smuce)

  ## Plot type 1 error for SMUCE vs TF
  w = 5;  h = 5.5   
  mar = c(4.5,4.5,2.3,0.5)
  pdf(file=file.path(outputdir,"type-1-error-smuce-by-n.pdf"),width=w, height=h)

  # Plot type 1 error profile  
    lty.t1 = 1
    pchs.t1 = c(16,17)
    lwd.t1 = 1 
    cex.t1 = 1
## par(mar=c(5.1,4.1,4.1,2.1))
    par(mar=mar)
    ylim=c(0,1)
    cols.t1 = c("red","black")  # brewer.pal(3,"Set1")
    col.lvl = 'lightgrey'
    lty.lvl = 1
    lwd.lvl = 2
   
    pchs=c(16,17,NA)
    cols=c('black','red', col.lvl)
    lwds = c(lwd.t1, lwd.t1, lwd.lvl)
    legends=c('TG','SMUCE','0.05 level')

    for(ii in 1:length(errs.list)){ 
      errs = errs.list[[ii]]
      # Make plot
      if(ii == 1){
          plot(errs ~ nlist, lty=lty.t1, pch=pchs.t1[ii],lwd=lwd.t1,type='o',
               cex=cex.t1,ylim=ylim,axes=F,xlab="n", col = cols.t1[ii],
               ylab = "Type I error")
        abline(h=0.05, col = col.lvl,lty=lty.lvl,lwd=lwd.lvl)
        title(main=bquote(Type~I~error))
#        text(x=300,y=0.07,labels="0.05")
        axis(1);axis(2)
        legend("topright", lty=c(1,1,1),pch=pchs,col=cols,legend=legends,lwd=lwds)
        ## legend("topright",lwd=lwd.lvl, lty=lty.lvl, col = col.lvl, legend = "0.05 level")
      } else {
        lines(errs ~ nlist, lty=lty.t1, pch=pchs.t1[ii], lwd=lwd.t1,type='o',cex=cex.t1, col = cols.t1[ii])
      }

      binom.std.err = unlist(Map(function(p,n) sqrt(p*(1-p)/n), errs, nlist))
      middle = errs
      upper = errs + 1.96*binom.std.err
      lower = errs - 1.96*binom.std.err
      errbar(x = nlist[length(nlist):1], 
           y = middle[length(middle):1], 
           yplus = upper[length(upper):1],
           yminus = lower[length(lower):1],
           pch=pchs.t1[ii],col=cols.t1[ii],errbar.col = cols.t1[ii],add=T)
    }
  graphics.off()
