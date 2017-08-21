# Make sure you're working from [dropboxfolder]/code
## source('funs.R')
## source('testfuns.R')
## source('dualPathSvd2.R')
## source('selectinf/selectiveInference/R/funs.inf.R')
## outputdir = "output"
## codedir = "."

outputdir = "~/Desktop"
library(genlasso)
library(RColorBrewer)

###########################################
### Generate data for declutter example ###
###########################################

# Form truth and noisy data (three linear segments)
  beta1 = seq(from=1,to=10,by=1)
  beta2 = -0.5*seq(from=11,to=20,by=1) + 15
  beta3 = seq(from=21,to=30,by=1) - 15
  beta0 = c(beta1,beta2,beta3)
  n = length(beta0)

# Example of linear trend filtering
  sigma = 1
  tf.order = 1
  consec = 2
                                        #set.seed(30)
  myseed = 79  # some good seeds: 53, 73, 69
  print(myseed)
  set.seed(myseed) #set.seed(53)

  y0 = beta0 + rnorm(length(beta0),0,sigma)
  D = makeDmat(n,type="tf",order=1)
  maxsteps = 20
  f0 = dualpathSvd2(y0, D, maxsteps = maxsteps)

  ## Collect Gammat at stop time
  bic   = get.modelinfo(f0,y0,sigma,maxsteps, stoprule = 'bic')$ic
  stop.time = which.rise(bic,consec=consec) - 1
  stop.time = pmin(stop.time,n-consec-1)

  if(!(stop.time+consec < maxsteps)){
    stop('bic rule hasnt stopped!')
  }

  Gobj.new.with.stoptime = getGammat.with.stoprule(obj=f0,y=y0,
                                     condition.step = stop.time+consec,
                                     stoprule = "bic", sigma=sigma, type='tf',
                                     consec=consec, maxsteps=maxsteps,D=D)
  G = Gobj.new.with.stoptime$Gammat
  u = Gobj.new.with.stoptime$u

  # Check correctness of polyhedron:
  polyhedron.checks.out = function(y0,G,u){
    all.nonneg = all((G%*%y0 >= u))
    if(all.nonneg){
      return(all.nonneg)
    } else {
      print("failed Gy >= u test")
      return(all.nonneg)
    }
  }
  #stopifnot(polyhedron.checks.out(y0,G.bic,u.bic))
  if(!polyhedron.checks.out(y0,G,u)){
      print("polyhedron is problematic")
      next
  }

  # Conduct tests and record pvals + stoptimes
  states = get.states(f0$action)

  # declutter the last states
  final.model = states[[stop.time]]
  final.model.cluttered = sort(final.model)
  final.model.decluttered = declutter(final.model)

  # test only the decluttered states, with /their/ adjusted contrasts
  final.models = list(final.model.cluttered,final.model.decluttered)

  Vs = Ps = Cs = list()
  for(kk in 1:2){
    final.model = final.models[[kk]]
    contrasts = list()
    pvals = coords = c()
    if(stop.time > 0){
      for(ii in 1:length(final.model)){
        this.sign = f0$pathobj$s[f0$pathobj$B == final.model[ii]]
        if(length(this.sign)==0){
          f1 = dualpathSvd2(y0, D, maxsteps = stop.time-1)
          this.sign = f1$pathobj$s[f1$pathobj$B == final.model[ii]]
        }


        v =     make.v.tf.fp(test.knot = final.model[ii],
                              adj.knot  = final.model,
                              test.knot.sign = this.sign,#f1$adj.knot,
                              D = D)
        contrasts[[ii]] = v

        mycoord = final.model[ii]
        mypval  = pval.fl1d(y0,G,dik=v,sigma,u=u)
        pvals[ii] = mypval
        coords[ii] = mycoord


      }
    }

    Vs[[kk]] = contrasts
    Ps[[kk]] = pvals
    Cs[[kk]] = coords
    if(kk==1)print("before")
    if(kk==2)print("after")
    names(pvals) = final.model
    print(pvals)
  }

####################################
## Plot two decluttering examples ##
####################################
lcol.final.model.knots = 'lightgrey'
lcol.test.knot = 'blue'
lcols.knot = c(lcol.final.model.knots,lcol.test.knot)
lwd.test.knot = 1.5
lwd.knot = 1
lwd.knots = c(lwd.knot,lwd.test.knot)
lty.knots = c(3,3)
lty.signal = lty.est = 1
pch.dat = 16
lcol.signal = 'red'
lcol.est = 'blue'
lwd.est = 2
lwd.signal=2
pcol.dat = "grey50"
## pcols.contrast = brewer.pal(n=3,name="Set2")
pcols.contrast.tf = brewer.pal(n=3,name="Set2")
pcols.contrast.tf[2] = "grey75"
## lwd.knots = c(1,1)
pchs.contrast.tf = c(17,15)
ylim=c(-5,20)
ylab = ""
xlab = "Location"
lcol.hline = "lightgrey"
xlim = c(0,length(y0))
filenames = c("tf-declutter-firstknot.pdf","tf-declutter-secondknot.pdf")
w=h= 5
mar=4*c(1,.5,.15,.5)
for(ii in 1:2){
for(kk in 2:1){
    if( (kk==1 & ii ==1) | (kk==1 &  ii==2) | (kk==2 &  ii==1) | (kk==2 &  ii==2)  ){
    v = (Vs[[kk]])[[ii]]
    pval = (Ps[[kk]])[ii]
    final.model = final.models[[kk]]

    if(kk==2) {
        pdf(file=file.path(outputdir, filenames[ii]),width=w,height=h)
        par(mar=mar)
        plot(NA, ylim = ylim, xlim=xlim, xlab=xlab,ylab=ylab,axes=F)
        axis(1);axis(2)
    }
    points(y0,pch=pch.dat, col = pcol.dat)
    abline(h=0,col=lcol.hline)
#      lines(beta0,col=lcol.signal, lwd=lwd.signal)
    points(v*7-kk*.3,col=pcols.contrast.tf[kk],pch=pchs.contrast.tf[kk]);
    abline(v=final.model[ii]+.5+.5, lwd=lwd.test.knot, col = lcol.test.knot,lty=lty.knots);
    abline(v=final.model+.5+.5, lwd=lwd.knot, col = lcol.final.model.knots,lty=lty.knots)
    lines(f0$beta[,stop.time+1], col = lcol.est, lwd = lwd.est)
    ## text(x=final.model[ii], y = -5-1.5*kk+1, label = paste(  (if(kk==1) "Original" else "Post-processed"),"p-value: ",signif(pval,3)))

        ## if(kk==2 & ii == 1){#if(kk==2 & ii == 1){
            legend("topleft", inset = 0,
                            lty = c(NA, lty.signal, lty.est, lty.knots[1],NA,NA)[-2],
                            lwd = c(NA, lwd.signal, lwd.est, lwd.knots[1],NA,NA)[-2],
                            pch = c(pch.dat, NA,NA,NA, pchs.contrast.tf, pchs.contrast.tf)[-2],
                            col = c(pcol.dat, lcol.signal,lcol.est,lcols.knot[1],pcols.contrast.tf)[-2],
                            legend=c("Data","Signal","Estimate","Knot","Original contrast", "Decluttered contrast")[-2],bg='white')
        ## }
        }
    }

    ## if(ii==2){
    ##     kk=2
    ##     v = (Vs[[kk]])[[ii]]
    ##     pval = (Ps[[kk]])[ii]
    ##     points(v*7-kk*.3,col=pcols.contrast.tf[kk],pch=pchs.contrast.tf[kk]);
    ## }
    graphics.off()
    }


#######################################################
### 1d Fused Lasso declutter (generate data + plot) ###
#######################################################

## Generate data and plot (in one swipe)
    set.seed(29)
    sigma = .5
    lev1=2; lev2=5; lev3=3; n = 60
    beta0 = rep(c(lev1,lev2,lev3),each=n/3)
    set.seed(0)

        y0    = beta0 + rnorm(n, 0, sigma)
        D = makeDmat(n,ord=0)
        maxsteps=10
        f0 = dualpathSvd2(y0,D,maxsteps=10)
        bic = get.modelinfo(f0,y0,sigma,D=D,maxsteps=10)$ic
        consec=2
        stop.time = which.rise(bic,consec)-1
        states = get.states(f0$action)

        final.model = states[[stop.time+1]]
        final.model.before.declutter = final.model
        final.model.after.declutter  = declutter(final.model)

    if(!(stop.time+consec < maxsteps)){
        stop('bic rule hasnt stopped!')
    }

    Gobj.new.with.stoptime = getGammat.with.stoprule(obj=f0,y=y0,
                                        condition.step = stop.time+consec,
                                        stoprule = "bic", sigma=sigma, type='tf',
                                        consec=consec, maxsteps=maxsteps,D=D)
    G = Gobj.new.with.stoptime$Gammat
    u = Gobj.new.with.stoptime$u

    # Check correctness of polyhedron:
    polyhedron.checks.out = function(y0,G,u){
        all.nonneg = all((G%*%y0 >= u))
        if(all.nonneg){
        return(all.nonneg)
        } else {
        print("failed Gy >= u test")
        return(all.nonneg)
        }
    }
    #stopifnot(polyhedron.checks.out(y0,G.bic,u.bic))
    if(!polyhedron.checks.out(y0,G,u)){
        print("polyhedron is problematic")
        next
    }

        # Two contrasts, before clustering
        v.before1.with.zeroes = getdvec(obj=f0,y=y0,k=2,klater=stop.time, type="segment")
        v.before1 = v.before1.with.zeroes
        v.before1[v.before1==0] = NA
        v.before2.with.zeroes = getdvec(obj=f0,y=y0,k=3,klater=stop.time, type="segment")
        v.before2 = v.before2.with.zeroes
        v.before2[v.before2==0] = NA

## p1.before = pval.fl1d(y0,G,dik=v.before1)
        p1.before = pval.fl1d(y0,G,dik=v.before1.with.zeroes,sigma=sigma,u=u)
        p2.before = pval.fl1d(y0,G,dik=v.before2.with.zeroes,sigma=sigma,u=u)

        # One contrast, after clustering
        f1 =  dualpathSvd2(y0,D,maxsteps=stop.time)
        this.sign = f1$pathobj$s[which(f1$pathobj$B == final.model.after.declutter[ii])]

        D = makeDmat(n,ord=0)
        v.after.with.zeroes = make.v.tf.fp(test.knot = final.model.after.declutter[ii],
                                adj.knot = final.model.after.declutter,
                               test.knot.sign = this.sign, D=D)
        v.after = v.after.with.zeroes
        vtol = 1E-10
        v.after[abs(v.after)<vtol] = NA
        p.after = pval.fl1d(y0,G,dik=v.after.with.zeroes,sigma=sigma,u=u)

## plot settings
w=5
h=5
pcols.contrast = brewer.pal(n=3,name="Set2")
pcols.contrast[3]=pcols.contrast[1]
pcols.contrast[2]='grey75'
ylim = c(-2,10)
xlim = c(0,60)
lty.knot = 3
lcol.knot = 'lightgrey'
lcol.test.knot = 'blue'
pchs.contrast = c(15,17,15)
lcol.est = 'blue'
lcol.hline = 'lightgrey'
filename = "1dfl-declutter-example.pdf"
mar=4*c(1,.5,.15,.5)
lwd.est = 2
lty.est=1

pdf(file.path(outputdir,filename), width=w,height=h)
par(mar=mar)
    plot(y0,pch=16, col = 'grey50',ylim=ylim,axes=F,xlab=xlab,ylab=ylab,xlim=xlim)
    axis(1);axis(2);
    lines(f0$beta[,stop.time+1],col=lcol.est,lwd=lwd.est,lty=lty.est)
    abline(v=states[[stop.time+1]],col=lcol.knot, lty = lty.knot)
    ## points(x=1:n, y=v.before1, col = pcols.contrast[1], pch=pchs.contrast[1])
    points(x=1:n, y=v.after*1.3,   col = pcols.contrast[3], pch=pchs.contrast[3])
    points(x=1:n, y=v.before2, col = pcols.contrast[2], pch=pchs.contrast[2])
    abline(h=0, col = lcol.hline)
    abline(v=final.model.before.declutter,col=lcol.knot, lty = lty.knot)
    abline(v=final.model.before.declutter[3], col = lcol.test.knot, lty = lty.knot,lwd=lwd.test.knot)

    ## text(x=36,y=-2, label = paste("Original p-value:",
    ##                             round(p1.before,3)))
    ## ## text(x=36,y=-2, label = paste("Original p-value:",
    ## ##                             round(p1.before,3),"(location 39)\n",
    ## ##                             "                                 ",
    ## ##                             round(p2.before,3), "(location 40)"))
    ## text(x=40,y=-2.7, label = paste("Post-processed p-value:", signif(p.after,3) ))
    legend("topleft",
        pch=c(pch.dat,NA,NA,pchs.contrast[c(2,3)]),
        lty=c(NA,lty.est,lty.knot,NA,NA),
        lwd=c(NA,lwd.est,lwd.knot,NA,NA),
        col = c(pcol.dat, lcol.est, lcol.knot, pcols.contrast[c(2,3)]),
        legend=c("Data", "Estimate","Knot",
                 "Original contrast", "Decluttered contrast"),
        bg="white")
    ## legend("topleft",
    ##     pch=c(pch.dat,NA,NA,pchs.contrast),
    ##     lty=c(NA,1,lty.knot,NA,NA,NA),
    ##     lwd=c(NA,2,lwd.knot,NA,NA,NA),
    ##     col = c(pcol.dat, lcol.est, lcol.knot, pcols.contrast),
    ##     legend=c("Data", "Estimate","knot","Original contrast (location 40)",
    ##                 "Original contrast (location 39)", "Contrast after declutter"),
    ##     bg="white")
graphics.off()
