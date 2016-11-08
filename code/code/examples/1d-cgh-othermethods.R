# working directory should be [genlassoinf]/code
workingdir = '~/../../media/shyun/Bridge/Dropbox/CMU/courses(CURRENT)/genlassoinf/code'
setwd(workingdir)
source("settings.R")
source('funs.R')
source('examples/testfuns.R')
source('dualPathSvd2.R')
lapply(c("genlasso","pryr"), require, character.only = TRUE)
library(changepoint)
library(wbs)
library(RColorBrewer)
#library(polypath)

#########################
# 2004 Hot spot paper ###
#########################

  library(cghFLasso)
  data(CGH)
  summary(CGH)
  w=h=5
  pdf(file=file.path(outputdir,"cgh-hotspotpaper2004.pdf"))
    CGH.FL.obj1<-cghFLasso(CGH$GBM.y)
    fit.2004 = CGH.FL.obj1[[1]]
    lcol.2004 = "red"
    lwd.2004 = 2
    dat.2004 = CGH.FL.obj1[[2]]
    pcol.2004 = "grey50"
    pch.2004 = 16
    plot(dat.2004, col = pcol.2004, pch = pch.2004,axes=F)
    axis(1);axis(2)
    lines(fit.2004, col = lcol.2004, lwd = lwd.2004)
    cp.2004 = which(sapply(1:(length(fit.2004)-1), function(ii) fit.2004[ii] != fit.2004[ii+1] ))
    abline(v=cp.2004, col = 'lightgrey',lty = 3)
  graphics.off()


############################
##### CP package ###########
############################

# Fit fused lasso paths
  y0 = CGH$GBM.y
  n = length(y0)

## PELT method
  m.pelt = cpt.mean(y0,method="PELT",class=F, penalty = "MBIC")
  m.pelt # plot(m.pelt)
  
## SBS method
  ss = sbs(y0)
  s.cpt <- unlist(changepoints(ss)$cpt.th)
  
## WBS
  w = wbs(y0)
  w.cpt = changepoints(w)$cpt.ic$ssic.penalty

## 1d sparse fused lasso (2004 paper)
  cp.2004
  
## Combine All
  cpts = list(pelt=m.pelt,sbs=s.cpt, w=w.cpt, sfl = cp.2004)
  cpts = lapply(cpts,sort)
  segs = function(y0,cpt){
    cpt = c(0,cpt,length(y0))
    inds = Map(function(a,b){(a+1):b}, cpt[1:(length(cpt)-1)], cpt[2:length(cpt)])
    means = lapply(inds, function(ind){rep(mean(y0[ind]) ,length(ind) )})
    means = unlist(means)
    return(means)
  }

  all.cpts = unique(unlist(cpts))
  common.cpts = all.cpts[sapply(all.cpts, function(this.cpt){ all( unlist (lapply(cpts,function(this.method.cpts) (this.cpt %in% this.method.cpts) ))) })]

## Plot all  
  pch=16
  titles = c("PELT","Standard Binary Segmentation","Wild Binary Segmentation","Fused Lasso (FDR)")
  methodnames = c("pelt","sbs","wbs","fl-fdr")
  xlab = ylab = ""
  cols=brewer.pal(9,"Set1")
  w=8; h=4;
  pcol='lightgrey'
  lcol.loc = cols[3]
  lcol.common.loc = cols[2]
  lcol.est = cols[1]
  lty.common.loc=1
  lty.loc=1
  for(jj in 1:4){
    pdf(file=file.path(outputdir,paste0("cgh-",methodnames[jj],".pdf")),width=w,height=h)
    par(mar=c(2.1,3.1,3.1,3.1))
    plot(y0,pch=pch,col=pcol,axes=F,ylab=ylab,xlab=xlab)
    axis(1); axis(2)
    title(main=titles[jj])
    mymean = segs(y0,cpts[[jj]])
    lines(mymean,lwd=2,col=lcol.est)
    this.method.cpts = cpts[[jj]]
    lcol.cpts = rep(lcol.loc,length(this.method.cpts))
    lcol.cpts[which(this.method.cpts %in% common.cpts)] = lcol.common.loc
    lty.cpts = rep(lty.loc,length(this.method.cpts))
    lty.cpts[which(this.method.cpts %in% common.cpts)] = lty.common.loc
    show.cpts = rep(F,length(this.method.cpts))
    show.cpts[which(this.method.cpts %in% common.cpts)] = T

    abline.vertical.segment =function(v,ylim,lty,col,lwd,show.cpts){
      for(ii in 1:length(v)){
        my.v = v[ii]
        if(ii %in% which(show.cpts)) next()
        lines(x = rep(my.v,2), y = ylim, col = col[ii], lwd=lwd, lty=lty[ii])
      }
    }
    abline.vertical.segment(v=this.method.cpts, col=lcol.cpts,lty=lty.cpts,ylim = c(-2.5,-2.2),lwd=1,show.cpts = show.cpts)
    abline.vertical.segment(v=this.method.cpts, col=lcol.cpts,lty=lty.cpts,ylim = c(-2.3,-2.0),lwd=1,show.cpts = !show.cpts)
    
    abline(h=0,col='grey50')
    if(jj==1){
      legend("topright",
            col = c(pcol,lcol.loc,lcol.common.loc,lcol.est),
            lty = c(NA,lty.loc,lty.common.loc,1),
            lwd = c(NA,1,1,2),
            pch = c(pch,NA,NA,NA),
            seg.len = 1,
            legend = c("Data","Locations","Common locations","Piece-wise means"),
            bg = "white")
    }
    graphics.off()
  }
  
#########################
## Sample Splitting #####
#########################
# 1st half: estimate the model
  y0 = CGH$GBM.y
  y1 = y0[2*(1:length(y0/2))]
  y2 = y0[2*(1:length(y0/2))-1]
  n1 = length(y1)
  maxsteps= 100
  D = rbind(getDmat(n/2,order=0))
  f0 = dualPathSvd2(y1,D,maxsteps=maxsteps)




