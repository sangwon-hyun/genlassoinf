## Make sure you're working from [dropboxfolder]/code
source("settings.R")
source('funs.R')
source('testfuns.R')
source('dualPathSvd2.R')
library(RColorBrewer)

outputdir = "output" # Ryan

########################################################
#### Generate p-values and simulation quantities #######
########################################################

  nsim = 500
  n = 60
  sigma = 1
  lev1= 0
  lev2= 2
 # lev2list=c(0,0.5,1,1.5,2)
  lev2list = c(0,1,2)*.5
  numsteps=1
#  spike   = introexample(testtype = "spike", nsim=nsim,sigma=sigma,lev1=lev1,lev2=lev2,lev2list=lev2list,numsteps=numsteps,verbose=T)
#  segment = introexample(testtype = "segment", nsim=nsim,sigma=sigma,lev1=lev1,lev2=lev2,lev2list=lev2list,numsteps=numsteps,verbose=T)
  #save(file = file.path(outputdir,"onejump-example-finersignal.Rdata"), list = c("lev1","lev2list","nsim","spike","segment","proportions0", "proportions1"))
  #save(file = file.path(outputdir,"onejump-example.Rdata"), list = c("lev1","lev2list","nsim","spike","segment","proportions0", "proportions1"))

###########################################
#### plotting QQ plots for p-values #######
###########################################
  load(file.path(outputdir,"onejump-example.Rdata"))

    set.seed(0)
    sigma = 1
    y0    = onejump.y(returnbeta=F,lev1=0,lev2=2,sigma=sigma,n=60)
    beta0 = lapply(c(0:2), function(lev2) onejump.y(returnbeta=T,lev1=0,lev2=lev2,sigma=sigma,n=60))
    beta0.middle = onejump.y(returnbeta=T,lev1=0,lev2=1,sigma=sigma,n=60)
    beta0.top = onejump.y(returnbeta=T,lev1=0,lev2=2,sigma=sigma,n=60)
    x.contrasts = c(1:60)
    v.spike = c(rep(NA,29),.5*c(-1,+1)-2, rep(NA,29))
    v.segment = c(.3*rep(-0.7,30)-2 , .3*rep(+0.7,30)-2 )

    xlab = "Location"
    w = 5; h = 5
    pch = 16; lwd = 2
    pcol = "gray50"
    ylim = c(-3,5)
    mar = c(4.5,4.5,0.5,0.5)
    xlim = c(0,70)

    xticks = c(0,2,4,6)*10
    let = c("A","B")
    ltys.sig = c(2,2,1)
    lwd.sig = 2
    pch.dat = 16
    pcol.dat = "grey50"
    pch.contrast = 17
    lty.contrast = 2
    lcol.sig = 'red'
    pcol.spike=3
    pcol.segment=4
    pcols.delta =   pcols.oneoff = brewer.pal(n=3,name="Set2")
    pch.spike = 15
    pch.segment = 17
    cex.contrast = 1.2

  ##################################
  ## Example of data and contrast ##
  ##################################
  pdf("output/onejump-example-data-and-contrast.pdf", width=5,height=5)

    par(mar=c(4.1,3.1,3.1,1.1))
      plot(y0, ylim = ylim,axes=F, xlim=xlim, xlab = xlab, ylab = "", pch=pch, col=pcol);
      axis(1, at = xticks, labels = xticks); axis(2)
      for(ii in 1:3) lines(beta0[[ii]],col="red",lty=ltys.sig[ii], lwd=lwd.sig)
      for(ii in 0:2) text(x=65,y=ii, label = bquote(delta==.(ii)))
      points(v.spike~x.contrasts, pch = pch.spike, col = 3)
      points(v.segment~x.contrasts, pch = pch.segment, col = 4)
      abline(h = mean(v.segment,na.rm=T), col = 'lightgrey')
      legend("topleft", pch=c(pch.dat,NA,pch.spike,pch.segment),
             lty=c(NA,1,NA,NA), lwd=c(NA,2,NA,NA),
             col = c(pcol.dat, lcol.sig, pcol.spike,pcol.segment),
             pt.cex = c(cex.contrast, NA, cex.contrast, cex.contrast),
             legend=c("Data", "Mean","Spike contrast", "Segment contrast"))
     title(main=expression("Data example"))
  graphics.off()

