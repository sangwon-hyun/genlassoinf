## Make sure you're working from [dropboxfolder]/code
source("settings.R")
source('funs.R')
source('testfuns.R')
source('dualPathSvd2.R')
library(genlasso)
library(RColorBrewer)
library(Hmisc)
outputdir = "output"

# Helper functions
  getpow = function(verdicts,igrain,loc){  return(sum(verdicts[igrain,,loc],na.rm=T)/pmax(1,sum(!is.na(verdicts[igrain,,loc]))))   }
  getpow.oracle = function(verdicts,igrain) { return(sum(verdicts[igrain,],na.rm=T)/pmax(1,sum(!is.na(verdicts[igrain,]))))  }

# Get powers
  load(file=file.path(outputdir, "paper-BIC-onejump-segmentsize10.Rdata"))
  loc=10
  load(file=file.path(outputdir, paste0("paper-BIC-onejump-oracle-segmentsize", 20/2, ".Rdata")))
  powers.10 = list(powers.bic = sapply(1:ngrain,function(igrain) getpow(verdicts.bic,igrain,loc)),
                powers.aic = sapply(1:ngrain,function(igrain) getpow(verdicts.aic,igrain,loc)),
                powers.fixed1 = sapply(1:ngrain,function(igrain) getpow(verdicts.fixed1,igrain,loc)),
                #powers.fixed2 = sapply(1:ngrain,function(igrain) getpow(verdicts.fixed2,igrain,loc)),
                powers.oracle = sapply(1:ngrain,function(igrain) getpow.oracle(verdicts.oracle,igrain)))
  verdicts.10 = verdicts.bic

  load(file=file.path(outputdir, "paper-BIC-onejump-segmentsize20.Rdata"))
  loc=20
  powers.20 = list(powers.bic = sapply(1:ngrain,function(igrain) getpow(verdicts.bic,igrain,loc)))
  verdicts.20 = verdicts.bic
 
  load(file=file.path(outputdir, "paper-BIC-onejump-segmentsize30.Rdata"))
  loc=30
  powers.30 = list(powers.bic = sapply(1:ngrain,function(igrain) getpow(verdicts.bic,igrain,loc)))
  verdicts.30 = verdicts.bic
 
  load(file=file.path(outputdir, "paper-BIC-onejump-segmentsize40.Rdata"))
  loc=40
  powers.40 = list(powers.bic = sapply(1:ngrain,function(igrain) getpow(verdicts.bic,igrain,loc)))
  verdicts.40 = verdicts.bic
  
##############  
# Plot data ##
##############  
  mar = c(4.5,4.5,2.5,0.5)
  w=h=5
  set.seed(0)
  pch.dat=16
  pcol.dat="grey50"
  lcol.sig="red"
  lwd.sig = 2
  lev1 = 0
  lev2 = 3
  n = 20
  sigma = 1
  ylab = ""
  xlab = "Coordinate"
  lty.sig = 1
  xlim = c(0,24)
  ylim = c(-2,5)

  pdf(file=file.path(outputdir,"bic-onejump-data.pdf"),width=w,height=h)
    par(mar=mar)
    beta0 = rep(c(lev1,lev2),each=n/2)
    y0    = beta0 + rnorm(n, 0,sigma)
    plot(y0, col = pcol.dat, axes=F,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch.dat)
    axis(1);axis(2)
    lines(beta0, col=lcol.sig,lwd= lwd.sig)
    legend("topleft",lty=c(lty.sig,NA), lwd=c(lwd.sig,NA), 
            col = c(lcol.sig,pcol.dat), pch = c(NA,pch.dat), 
            legend = c("Signal","Data"))
    title(main="Example of data and signal \n(n=20)", line=1)
    text(x=22,y=3.2,labels=bquote(sigma==.(sigma)))
    text(x=22,y=2.8,labels=bquote(delta==.(lev2)))
    abline(v=10,col='lightgrey',lty=3)
  graphics.off()    

################
# Plot powers ##
################  
  nn=4
  cols = brewer.pal(nn,"Set1")
  lwd=2
  pchs = c(1:nn)
  ylab = "Power"
  ylim=c(0,1)

  mar = c(4.5,4.5,3.5,0.5)
  w=5;h=5.5


  #####################
  #### within n=20 ####
  #####################
  pdf(file=file.path(outputdir,"bic-onejump-powers.pdf"),width=w,height=h)
    par(mar=mar)
    matplot(y=do.call(cbind,powers.10),x=lev2list,lty=1,col=cols,lwd=lwd,type='o',
        pch=pchs,ylab=ylab, xlab=bquote(delta),axes=F,ylim=ylim)
    axis(1);   axis(2)
    legend("bottomright",lwd=rep(lwd,nn), lty=rep(1,nn), pch = pchs[c(2,1,3,4)] ,col=cols[c(2,1,3,4)], legend = c("AIC","BIC","Fixed","Oracle")) 
    title(main=bquote(atop(Power,(location==10))),line=1.5)
    #title(main=bquote(atop(Conditional~Powers,(location==10~~with~n==20))),line=1.5)

    # binomial confidence bands
    for(method.i in 1:length(powers.10)){
      pp = powers.10[[method.i]]
      binom.band.width=rep(NA,ngrain)
      ns = rep(NA,ngrain)
      for(igrain in 1:ngrain){
        n.effective = sum(apply(verdicts.10[igrain,,],1,function(myrow)sum(!(is.na(myrow))))!=0)
        print(n.effective)
        ns[igrain] = n.effective
        binom.band.width[igrain] = 1.96*sqrt(pp[igrain] * (1-pp[igrain])/n.effective)
      }
      errbar(x = lev2list,
             y = pp,
             yplus = pp + binom.band.width,
             yminus = pp - binom.band.width, add=T, col=cols[method.i],errbar.col = cols[method.i],pch = pchs[method.i])
    }
    
  graphics.off()
  
  

  ###########################
  ### Among n=20,40,60,80 ###
  ###########################
  lwd.n20 = 3
  cex=.5
  ylim = c(0,1)
  cols.by.n = c(1,2,1,2)
  mar = c(4.5,4.5,3.5,0.5)
  w=5;h=5.5
pdf(file=file.path(outputdir,"bic-onejump-powers-by-n.pdf"),width=w,height=h)
par(mar=mar)
    plot(y=powers.10$powers.bic,x=lev2list,lty=1,col=cols.by.n[1],lwd=lwd.n20,type='o',pch=pchs[1],ylab=ylab, xlab=bquote(delta),axes=F,ylim=ylim,cex=cex)
    lines(y=unlist(powers.20),x=lev2list, col = cols.by.n[2], lwd = lwd, pch=pchs[1], type='o', lty=2,cex=cex)
    lines(y=unlist(powers.30),x=lev2list, col = cols.by.n[3], lwd = lwd, pch=pchs[1], type='o', lty=3,cex=cex)
    lines(y=unlist(powers.40),x=lev2list, col = cols.by.n[4], lwd = lwd, pch=pchs[1], type='o', lty=4,cex=cex)
    axis(1);   axis(2)
    ## title(main="Conditional Powers\n Testing location n/2") #
    title(main=bquote(atop(Power,(location==n/2))), line=1.5)
    legend("bottomright",lwd=rev(c(lwd.n20,rep(lwd,nn-1))), lty=rev(seq(1,4)), col=rev(cols.by.n), legend = rev(c("BIC (n=20)",
                                                                                       "BIC (n=40)",
                                                                                       "BIC (n=60)",
                                                                                       "BIC (n=80)")))
                                                                                       
    powers.bic = list(powers.10$powers.bic, unlist(powers.20),unlist(powers.30),unlist(powers.40))    
    verdicts.list = list(verdicts.10,verdicts.20,verdicts.30,verdicts.40)
    # binomial confidence bands
    for(size.i in 1:length(powers.bic)){
      pp       = powers.bic[[size.i]]
      verdicts = verdicts.list[[size.i]]
#      for(power.i in 1:length(powers.bic)){
      binom.band.width=rep(NA,ngrain)
      ns = rep(NA,ngrain)
      for(igrain in 1:ngrain){
        n.effective = sum(apply(verdicts[igrain,,],1,function(myrow)sum(!(is.na(myrow))))!=0)
        ns[igrain] = n.effective
        binom.band.width[igrain] = 1.96* sqrt(pp[igrain] * (1-pp[igrain])/n.effective)
      }
      errbar(x = lev2list,
             y = pp,
             yplus = pp + binom.band.width,
             yminus = pp - binom.band.width, add=T, col=cols.by.n[size.i],errbar.col = cols.by.n[size.i],pch = NA)
#    }
    }
                                                                                       
  graphics.off()



############################
## See average stoptimes ###
############################
for(stoptimes in list(stoptimes.aic, stoptimes.bic)){
  print(round(apply(stoptimes,1,function(myvec)c(mean(myvec,na.rm=T),sd(myvec,na.rm=T))),2))
}
round(lev2list,2)




## Aggregate everything for n=20 (slower but readable)
#  powers.20 = list(powers.bic  = getpowers.from.chunks(20, "verdicts.bic"),
#                   powers.bic.decluttered  = getpowers.from.chunks(20, "verdicts.bic.decluttered")
#                   powers.ebic = getpowers.from.chunks(20, "verdicts.ebic"),
#                   powers.aic = getpowers.from.chunks(20, "verdicts.aic"),
#                   powers.fixed1 = getpowers.from.chunks(20, "verdicts.fixed1"),
#                   powers.fixed2 = getpowers.from.chunks(20, "verdicts.fixed2"),
#                   powers.oracle = getpowers.oracle(20))
#                   
#  stoptimes.20 = getstoptimes(n=20,ngrain=20,nsim=100000,nchunks=100,file.id="bic-onejump-segmentsize-allbics")

# # save(powers.20, stoptimes.20, ngrain,nsim,nchunks,sigmalist, file=file.path(outputdir, "bic-info-20.Rdata"))
#  load(file=file.path(outputdir, "powers-20.Rdata"))

### make plot of powers
##  pdf(file.path(outputdir,paste0( file.id,"-largesim", n/2, ".pdf")), width=10, height=6)
#  pdf(file.path(outputdir,paste0( file.id,"-feb2016", n/2, ".pdf")), width=10, height=6)
#    makeplot(20, sigmalist, powers.20, stoptimes.20)
#  dev.off()


#  




### Seeing why BIC2 beats BIC
#  n = 20
#  load(file=file.path(outputdir, paste0("bic-onejump-segmentsize-allbics", n/2, ".Rdata")))
#  inds = inds.1
#  noisei = 2
#  step=2
#  table(round(ssl.bic2[[noisei+1]]))
#  inds = which(stoptimes.bic2.cond[noisei,]==step+1)
#  pvalslist = list()
#  for(jj in 1:length(inds)){
#    tempmat = rbind(signif(pvals.bic2[noisei,inds[jj],],2), signif(pvals.fixed1[noisei,inds[jj],],2))
#    rownames(tempmat) = c(stoptimes.bic2.cond[noisei,inds[jj]], 1)
#    stoptimes.bic2.cond[,inds[jj]]
#    pvalslist[[jj]] = tempmat
#  }




## overlay JUST the BIC curves (and the two oracles)
#  pdf(file.path(outputdir, "bic-overlayed-onejump.pdf"), width = 8, height=8)
#  nlist = 2*c(10,40,90) #n = 10  
#  for(zz in 1:3){
#    zzold2 = zz
#    n = nlist[zz]
#    cat('\n', n, "out of", nlist)
#    source("settings.R")
#    load(file=file.path(outputdir, paste0("bic-onejump-segmentsize", n/2, ".Rdata")))
#    source("settings.R")
#    zz = zzold2
#    rm(zzold2)
#    print(zz)
#  # plot powers
#    loc = n/2
#    xlim = c(0,max(sigmalist))
#    ylim = c(0,1)
#    plot(powers.bic[,loc] ~ sigmalist, type = 'l', lwd=2, lty = zz,xlim=xlim, ylim=ylim,axes=F)
#    title("Conditional Power at correct location (n/2), for one-jump")
#    axis(1, padj = zz); axis(2);

#  # plot oracle
#    lines(powers.oracle~sigmalist, col = 'blue', lwd=2, lty = zz)
#    par(new=T)
#  }
#    legend("topright", lty = c(1:3), legend = nlist/2)
#  dev.off()
#  
#  
#  

## overlay JUST the BIC curves on SAME SCALE (and the two oracles)
#  pdf(file.path(outputdir, "bic-overlayed-onejump-sameaxis.pdf"), width = 8, height=8)
#  nlist = 2*c(10,40,90) #n = 10  
#  xlim = c(0,6)
#  plot(NA, type = 'l', lwd=2, ylim = c(0,1), xlim = xlim, axes = F, xlab=expression(sigma), ylab = "condit. powers", lty = zz)
#  title("Conditional Power at correct location (n/2), for one-jump")
#  axis(1); axis(2);
#  for(zz in 1:3){
#    n = nlist[zz]
#    cat('\n', n, "out of", nlist)
#    source("settings.R")
#    zzold2 = zz
#    load(file=file.path(outputdir, paste0("bic-onejump-segmentsize", n/2, ".Rdata")))
#    source("settings.R")
#    zz = zzold2
#    rm(zzold2)
#  # plot powers
#    loc = n/2
#    lines(powers.bic[,loc] ~ sigmalist, type = 'l', lwd=2, lty = zz)

#  # plot oracle
#    lines(powers.oracle~sigmalist, col = 'blue', lwd=2, lty = zz)
#  }
#  legend("topright", lty = c(1:3), legend = nlist/2)
#  dev.off()
#  
#  
#  
## plot data example
#  load(file=file.path(outputdir, "bic-onejump.Rdata"))
#  source("settings.R")
#  pdf(file.path(outputdir,"bic-onejump-example.pdf"), width=12, height=5)
#  par(mfrow = c(2,3))
#    for(n in c(20,80)){
#    # plot data example
#    sigma = .1
#    beta0 = onejump.y(returnbeta=T, lev1=0, lev2=2, sigma=sigma, n=n)  # this could change
#    y0    = onejump.y(returnbeta=F, lev1=0, lev2=2, sigma=sigma, n=n)
#    plot(y0,xlab="",ylab="");lines(beta0,col='red');title(paste("example data with noise=",sigma))
#    
#    sigma = .5
#    beta0 = onejump.y(returnbeta=T, lev1=0, lev2=2, sigma=sigma, n=n)  # this could change
#    y0    = onejump.y(returnbeta=F, lev1=0, lev2=2, sigma=sigma, n=n)
#    plot(y0,xlab="",ylab="");lines(beta0,col='red');title(paste("example data with noise=",sigma))
#    
#    sigma = 1
#    beta0 = onejump.y(returnbeta=T, lev1=0, lev2=2, sigma=sigma, n=n)  # this could change
#    y0    = onejump.y(returnbeta=F, lev1=0, lev2=2, sigma=sigma, n=n)
#    plot(y0,xlab="",ylab="");lines(beta0,col='red');title(paste("example data with noise=",sigma))
#    }
#  dev.off()
#  
#  
#  

## Seeing binomial confidence band; it is very small.
#  n=20
#  load(file=file.path(outputdir, paste0("bic-onejump-segmentsize", n/2, ".Rdata")))
#  ls()
#  sigmalist[12]
#  sd(verdicts.bic.naive[12,,10],na.rm=T)/100
#  mean(verdicts.bic.naive[12,,10],na.rm=T)
#  mean(verdicts.proxaic[12,10],na.rm=T)
#  powers.proxaic[12]
#  ind.ends = seq(from=1,to=10000,by=100)
#  
## Find out why there isn't a boost in power.  
#  n=20
#  load(file=file.path(outputdir, paste0("bic-onejump-segmentsize", n/2, ".Rdata")))
#  ii=5
#  for(isim in 1:100){
#    loc = n/2
#    proxwidth = round(log(n)) #2
#    proxlocs = (loc-proxwidth):(loc+proxwidth)
#    print(verdicts.bic[ii,isim,proxlocs])
#      readline() 
#    print(verdicts.bic[ii,isim,loc])
#      readline() 
#  }
#  
## find out why there isn't a boost in power; take 2
#  for(n in 2*c(10,40,90)[1]){
#  load(file=file.path(outputdir, paste0("bic-onejump-segmentsize-allbics", n/2, ".Rdata")))
#  ii=3
#  loc = n/2
#  proxwidth = n*0.15#round(log(n)) #2
#  cat("width is ", proxwidth, fill=T)
#  proxlocs = (loc-proxwidth):(loc+proxwidth)
#  verdicts1 = c()
#  verdicts2 = verdicts.ebic[ii,,loc]
#  for(isim in 1:nsim){
#    verdicts = verdicts.ebic[ii,isim,proxlocs]
#    verdicts1[isim] =  (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
#  }
#  
#  readline("proximate version")
#  num1 = sum(verdicts1,na.rm=T)
#  denom1 = sum(!is.na(verdicts1))
#  print(num1)
#  print(denom1)
#  print(num1/denom1)

#  readline("nonproximate version")
#  num2 = sum(verdicts2,na.rm=T)
#  denom2 = sum(!is.na(verdicts2))
#  print(num2)
#  print(denom2)
#  print(num2/denom2)
#  }
#  sigmalist[ii]
#  
## why does proxbic not do better?
###?



## Add back the null cases to the EBIC criterion; to see what is swept under the rug by EBIC condit power ignoring the (increased) null cases compared to BIC

### calculate condit power at each correct jump location
#  nlist = 2*c(10,40,90) #n = 10  
#  for(n in nlist[1]){
#    cat('\n', n, "out of", nlist)
##    source("settings.R")
#    load(file=file.path(outputdir, paste0("bic-onejump-segmentsize-allbics", n/2, ".Rdata")))
# #   source("settings.R")
#    powers.bic = powers.bic.naive = 
#    powers.ebic = powers.ebic.naive = 
#    powers.oracle = powers.proxbic = powers.proxbic2 = powers.proxebic = powers.proxsbic = powers.proxaic = array(NA, ngrain)
#    # For each location,
#    for(ii in 1:ngrain){
#      # get exact condit powers
#      getpow = function(verdicts,ii,loc){  return(sum(verdicts[ii,,loc],na.rm=T)/pmax(1,sum(!is.na(verdicts[ii,,loc]))))   }
#      getpow.augmented = function(verdicts,ii,loc){  return(sum(verdicts[ii,,loc],na.rm=T)/pmax(1,sum(!is.na(verdicts[ii,,loc]))))   }
#      for(loc in 1:n){
#        powers.bic[ii,loc]        = getpow(verdicts.bic,ii,loc)
#        powers.ebic[ii,loc]       = getpow(verdicts.ebic,ii,loc)
#        powers.ebic.augmented[ii,loc] = getpow(verdicts.ebic.naive,ii,loc)
#      }
#      
#      # get approximate powers at true break coordinate
#      loc = n/2
#      proxwidth = .15*n#log(n) #2
#      proxlocs = (loc-proxwidth):(loc+proxwidth)
#      proxbic.verdict = proxbic2.verdict = proxebic.verdict = proxsbic.verdict = proxaic.verdict = c()
#      for(isim in 1:nsim){
#       # bic
#        verdicts = verdicts.bic[ii,isim,proxlocs]
#        proxbic.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
#       # bic2
#        verdicts = verdicts.bic2[ii,isim,proxlocs]
#        proxbic2.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
#       # ebic
#        verdicts = verdicts.ebic[ii,isim,proxlocs]
#        proxebic.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
#       # sbic
#        verdicts = verdicts.sbic[ii,isim,proxlocs]
#        proxsbic.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
#       # aic
#        verdicts = verdicts.aic[ii,isim,proxlocs]
#        proxaic.verdict[isim] = (if(!all(is.na(verdicts))) any(verdicts,na.rm=T) else NA)
#      }
#      powers.proxbic[ii] = sum(proxbic.verdict, na.rm=T)/pmax(sum(!is.na(proxbic.verdict)),1)
#      powers.proxbic2[ii] = sum(proxbic2.verdict, na.rm=T)/pmax(sum(!is.na(proxbic2.verdict)),1)
#      powers.proxebic[ii] = sum(proxebic.verdict, na.rm=T)/pmax(sum(!is.na(proxebic.verdict)),1)
#      powers.proxsbic[ii] = sum(proxsbic.verdict, na.rm=T)/pmax(sum(!is.na(proxsbic.verdict)),1)
#      powers.proxaic[ii] = sum(proxaic.verdict, na.rm=T)/pmax(sum(!is.na(proxbic.verdict)),1)

#      # oracle power at true break coordinate
#      powers.oracle[ii] = sum(verdicts.oracle[ii,],na.rm=T)/nsim
#    }
#    obj.list2 = c("powers.bic", "powers.bic.naive", 
#                  "powers.bic2", "powers.bic2.naive",
#                  "powers.ebic", "powers.ebic.naive",
#                  "powers.sbic", "powers.sbic.naive",
#                  "powers.aic", "powers.aic.naive",
#                  "powers.fixed1", "powers.fixed2", "powers.oracle",
#                  "powers.proxbic", "powers.proxbic2","powers.proxaic")
#    save(list=c(obj.list1,obj.list2), file=file.path(outputdir, paste0("bic-onejump-segmentsize-allbics", n/2, ".Rdata")))
#  }




### Some simpler plotting code
##plot(powers.bic.20$powers[,10], type = 'l')
##lines(powers.bic2.20$powers[,10],col='blue', type = 'l')
##plot(powers.bic.80$powers[,40], type = 'l')
##lines(powers.bic2.80$powers[,40],col='blue', type = 'l')
### Single location and prox tests have essentially identical power
##plot(powers.bic.20$powers[,10],type='l')
##lines(powers.bic.20$powers.prox,col='red')

