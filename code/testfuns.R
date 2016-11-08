## functions used for test*.R files; mostly generating toy data and producing p-values
onejump.y <- function(returnbeta=F, lev1=1, lev2=5, sigma=.5, n=60){
  beta = rep(c(lev1,lev2),each=n/2)
  y <- beta + rnorm(n,sd=sigma)
  return(if(returnbeta) beta else y)
}

twojump.y <- function(returnbeta=F,  n = 60,  sigma = .5){
  beta = rep(c(1,6,8),each=n/3)
  y <- beta + rnorm(n,sd=sigma)
  return(if(returnbeta) beta else y)
}

fourjump.y <- function(returnbeta = F,  n = 100,  sigma = 1,levs){
  if(length(levs)!=5) stop("Provide five levels!")
  beta = rep(levs,each=n/5)
  set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
  y = beta + rnorm(n,sd=sigma)
  return(if(returnbeta) beta else y)
}

alternjump.y <- function(returnbeta = F, lev1 = 1, lev2 = 5,  n = 60,  sigma = .5, seed = NA){
  beta = rep(c(lev1,lev2,lev1),each=n/3)
  if(is.na(seed)) {set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
  } else {set.seed(seed)}
  y = beta + rnorm(n,sd=sigma)
  return(if(returnbeta) beta else y)
}


#fourjump.y <- function(returnbeta = F, lev1 = 1, lev2 = 5,  n = 60,  sigma = .5, seed = NA){
#  beta = c(rep(2,n/5),rep(16,n/5), rep(8,n/5),rep(10,n/5),rep(4,n/5))
#  if(is.na(seed)) {set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
#  } else {set.seed(seed)}
#  y = beta + rnorm(n,sd=sigma)
#  return(if(returnbeta) beta else y)
#}



#  pdf("onejump.pdf")
#    plot(onejump.y(returnbeta=F))
#    lines(onejump.y(returnbeta=T), col = 'red')
#  dev.off()

#  pdf("twojump.pdf")
#    plot(twojump.y(returnbeta=F), ylab = "y", xlab = "x")
#    lines(twojump.y(returnbeta=T), col = 'red')
#  dev.off()

#  pdf("fourjump.pdf")
#    plot(fourjump.y(returnbeta=F), xlab = "x", ylab = "y")
#    lines(fourjump.y(returnbeta=T), col = 'red')
#  dev.off()

#  pdf("alternjump.pdf")
#    plot(alternjump.y(returnbeta=F), xlab = "x", ylab = "y")
#    lines(alternjump.y(returnbeta=T), col = 'red')
#  dev.off()

## collect pvalue and jump along the way
getpj <- function(k=1, geny, lev1, lev2,  sigma = .5){
  y <- geny(lev1=lev1,lev2=lev2)
  pval = jump = c()
  for(kk in 1:k){
    f1 <- fl1d(y,kk)
    G <- getGammat(f1,y,kk) # Gamma matrix after kk'th step
    d <- getdvec(f1,y,kk)   # dik vector after kk'th jump
    pval[kk] <- pval.fl1d(y,G,d,sigma)  
    jump[kk] <- as.integer(f1$B[kk])
  }
  return(c(pval, jump)) # pvalue of the real jump  
}

## generate null p's
## k0 is *location* of gap to test, k2 steps taken in algorithm
nullp <- function(k0 = 1, k2 = 5,  n = 60,  sigma = .5, type = "spike",approx=T, approxtype = "gsell", threshold=T){ 
  if( k0 <= 0 || k0 >= n || k2 <= 0 || k2 >= n ) stop("k0 and k2 out of bounds!")
  beta0 = rep(5,n)
  y <- beta0 + rnorm(n,sd=sigma)
  f <- fl1d(y,k2)
  G <- getGammat(f,y,k2)
  d1_spike <- rep(0,n); d1_spike[k0] <- -1; d1_spike[k0+1] <- +1; 
  d1_segment <- rep(0,n)
  # assume you're breaking the whole theta
  d1_segment[1:k0] <- (n - k0)/n
  d1_segment[(k0+1):n] <- -k0/n
  sk <- 1 # any sign is fine for this purpose
  d1_segment <- -sk *d1_segment
  
  if(type == "spike"){
    return(pval.fl1d(y, G, d1_spike, sigma, approx=approx, approxtype = approxtype, threshold=threshold))
  } else if (type == "segment"){
    return(pval.fl1d(y, G, d1_segment, sigma, approx=approx, approxtype = approxtype, threshold=threshold))
  } else {
    stop("Not coded yet")
  }
}

## generate non-null p's
## k0 is *location* of gap to test, k2 steps taken in algorithm, geny generates data.
nonnullp <- function(k0 = 30, k2 = 5, geny = onejump.y){ 
  y <- geny()
  f <- fl1d(y,k2)
  G <- getGammat(f,y,k2)
  d <- rep(0,n); d[k0] <- -1; d[k0+1] <- +1;
  return(pval.fl1d(y, G, d, sigma))
}

## record p-values for 1-jump situation
record = function(numsteps=5, lev1 = 1, lev2 = 5,  sigma = .5){
  y0 <- onejump.y(lev1 = lev1, lev2 = lev2)
  D <- dual1d_Dmat(length(y0))
  k0 <- 30
  f0 <- fl1d(y0,maxsteps=numsteps)#dualpathSvd2(y0,D,maxsteps=5)#fl1d(y0,maxsteps=5)
  G  <- getGammat(f0,y0,numsteps)
  pval = jump = c()
  for(kk in (-5):5){
    d1 <- rep(0,length(y0)); d1[k0+kk] <- -1; d1[k0+kk+1] <- +1;
    pval[kk+6] <- pval.fl1d(y0,G,d1,sigma)
  }
  names(pval) = (-5):5
  return(pval)
}


# plot results from test-hybridsegment.
# |filename|: name of Rdata file to read in
# |plotname|: name of pdf file to output
# ylab: what should the ylabel be?
plotpower = function(filename,plotname,ylab, width=10,height=6){

  load(file=file.path(outputdir,filename))
  pdf(file.path(outputdir,plotname),width=width,height=height)
   
    par(mfrow = c(1,2))
    # plot results from conditioning on step 1
    plot(x = levs, y = madecut[,1], ylim = c(0,1), ylab = ylab, col = "black", type = 'l')
    lines(x = levs, y = madecutsegment[,1], col = "red", lty=1)
    lines(x = levs, y = madecuthybridsegment[,1], col = 'blue', lty=1)
 
    # plot results from conditioning on step 2
    lines(x = levs, y = madecut[,2], col = 'black', lty = 2)
    lines(x = levs, y = madecutsegment[,2], col = 'red', lty = 2)
    lines(x = levs, y = madecuthybridsegment[,2], col = 'blue', lty=2)

    # make legend
    legend("topleft", col = rep(c("black","red","blue"),each=2), lty = c(1,2,1,2,1,2), legend = c("1spike","2-1spike","1segment","2-1segment","1hybridsegment","2-1hybridsegment"))

    # calculate proportions
    prop.spike.correct1 =  sapply(info, function(lst)sum(lst$correct1))/nsim # ryan told me not to do this.
    prop.segment.correct1 = sapply(info, function(lst)sum(lst$segmentcorrect1))/nsim
    prop.hybridsegment.correct1 = sapply(info, function(lst)sum(lst$hybridsegmentcorrect1))/nsim

    prop.spike.correct21 =  sapply(info, function(lst)sum(lst$correct21))/nsim # ryan told me not to do this.
    prop.segment.correct21 = sapply(info, function(lst)sum(lst$segmentcorrect21))/nsim
    prop.hybridsegment.correct21 = sapply(info, function(lst)sum(lst$hybridsegmentcorrect21))/nsim

    # plot proportions
    widths = levs-lev1
    plot(x = widths, y = prop.spike.correct1, ylim = c(0,1), ylab = "proportion of false nulls", xlab = "humpwidth", main = "proportion of false nulls by width of hump, \n at first step", type = 'l', cex=.5)
    lines(x = widths, y = prop.segment.correct1, col = 'red', type='l', cex=.5)
    lines(x = widths, y = prop.hybridsegment.correct1, col = 'blue', type='l', cex=.5)
    
    lines(x = widths, y = prop.spike.correct21, ylim = c(0,1), ylab = "proportion of false nulls", xlab = "humpwidth", main = "proportion of false nulls by width of hump, \n at first step", type = 'l', cex=.5, lty=2)
    lines(x = widths, y = prop.segment.correct21, col = 'red', type='l', cex=.5, lty=2)
    lines(x = widths, y = prop.hybridsegment.correct21, col = 'blue', type='l', cex=.5, lty=2)
    
    
    # make legend
    legend("bottomright",col = c("black","red","blue","black","red","blue"),lty = c(1,1,1,2,2,2), legend = c("spike1","segment1", "hybridsegment1","spike21","segment21", "hybridsegment21"))

  dev.off()
  }
  

# Function to collect pvals, for an introductory example in Example section
introexample = function(testtype = c("spike","segment"), nsim=nsim,sigma=sigma,lev1=lev1,lev2=lev2,lev2list=lev2list, numsteps=numsteps,verbose=T){

  testtype <- match.arg(testtype)

  pvals.correctlist = list()
  pvals.oneoff = c()

  # collect p values when detection is correct
  if(verbose) cat('\n','collecting p-values for correct location', '\n')
  for(ii in 1:length(lev2list)){
    templev2 = lev1+lev2list[ii]
    cat("\n","lev1 is ", lev1, "temporary lev2 is ", templev2, "\n")
    alldone = F
    jj = 0
    pvals.correct = c()
    while(!alldone){
      y0     = onejump.y(returnbeta=F,lev1=lev1,lev2=templev2,sigma=sigma,n=n)
      path   = dualpathSvd2(y0,dual1d_Dmat(n),maxsteps=numsteps,approx=T)
      G      = path$Gammat[1:path$nk[1],]
      d      = getdvec(obj=path, y=y0, k=1, usage = "dualpathSvd", type=testtype, matchstep=TRUE)
      if(path$pathobj$B[1] == n/2){
        pvals.correct[jj] = pval.fl1d(y=y0, G=G, dik=d, sigma=sigma, approx=TRUE, threshold=TRUE, approxtype="rob")
        jj = jj + 1
        if(verbose) cat("\r", jj, "of", nsim)
      }
      alldone = (jj == nsim)
    }
    pvals.correctlist[[ii]] = pvals.correct
  }
  
  if(verbose) cat('\n')
  
  # collect p values from when detection is one off
  if(verbose) cat('\n','collecting p-values for 1 off','\n')
  alldone=F
  jj=0
  ii=0
  kk=0
  while(!alldone){
    y0     = onejump.y(returnbeta=F,lev1=lev1,lev2=lev2,sigma=sigma,n=n)
    path   = dualpathSvd2(y0,dual1d_Dmat(n),maxsteps=numsteps,approx=T)
    G      = path$Gammat[1:path$nk[1],]
    d      = getdvec(obj=path, y=y0, k=1, usage = "dualpathSvd", type=testtype, matchstep=TRUE)
    if(abs(path$pathobj$B[1] - n/2) == 0){ ii = ii+1 }
    if(abs(path$pathobj$B[1] - n/2) == 1){
      pvals.oneoff[jj] = pval.fl1d(y=y0, G=G, dik=d, sigma=sigma, approx=TRUE, threshold=TRUE, approxtype="rob")
      jj = jj+1
      if(verbose) cat("\r", jj, "of", nsim)
    }
    kk = kk+1
    alldone = (jj == nsim)
  }
  
  if(verbose) cat('\n')
  
  return(list(pvals.correctlist = pvals.correctlist, 
              pvals.oneoff  = pvals.oneoff, 
              lev1 = lev1,
              lev2 = lev2,
              lev2list = lev2list,
              nsim = nsim,
              sigma = sigma,
              kk = kk,
              jj = jj,
              prop.oneoff = kk/jj))
}
  

# getting the proportions of correct detections, one off detections, and more off. (now obsolete)
getprops = function(lev1, lev2, nsim = 1000, n,sigma){
  firstvar = c()  
   for(jj in 1:nsim){
    y0     = onejump.y(returnbeta=F,lev1=lev1,lev2=lev2,sigma=sigma,n=n)
    path   = dualpathSvd2(y0,dual1d_Dmat(n),maxsteps=numsteps,approx=T)
    firstvar[jj] = path$pathobj$B[1]
  }
  firstvar = as.factor(firstvar)
  result = c(length(firstvar[firstvar==30])/nsim,
             length(firstvar[firstvar%in%c(29,31)])/nsim,
             length(firstvar[firstvar%in%c(26:28,32:34)])/nsim)
  names(result) = c("correct", "oneoff", "moreoff")
  return(result)
}






# Function that simulates brownian motion trajectory
# Taken from http://www.r-bloggers.com/generate-stock-option-prices-how-to-simulate-a-brownian-motion/
genstock = function(size, myseed=NA){
 if(!is.na(myseed)){set.seed(myseed)}
  y = rnorm(size)
  x = y
  for (i in 1:size){ x[i] = 1/sqrt(size)*(sum(y[1:i])*sqrt(i)) }
  return(x)
}




  
  
  
  
## Forms contrast using LRT  
  make.v.tf = function(test.knot, adj.knot, sign.test.knot, n, tf.order, approx=F){
    H = getH.trendfilter(n,1)              
    test.coord = test.knot+1
    adj.coord = c(1:(tf.order+1), (1+adj.knot))
    adj.coord = adj.coord[adj.coord!=test.coord]
    if(approx) adj.knot = c(max(adj.coord[adj.coord<test.coord]), 
                            min(adj.coord[adj.coord>test.coord]))
    adj.knot = adj.knot[abs(adj.knot) != Inf]
    X = H[, adj.coord]
    PH = fitted(lm(diag(n) ~ X -1))
    orthPH = (diag(n) - PH)
    v = as.numeric(H[,test.coord] %*% orthPH)
    v = v * sign.test.knot
    return(v)
  }

## Forms contrast using the best fit continuous line in two adjacent segments
  make.v.tf.direct = function(test.knot, adj.knot, sign.test.knot, n, tf.order){
    adjacent.knots = suppressWarnings(c(max(adj.knot[adj.knot<test.knot]), 
                                        min(adj.knot[adj.knot>test.knot])))
    left.knot = adjacent.knots[1]
    right.knot = adjacent.knots[2]
    if(right.knot == Inf) right.knot = n
    if(left.knot == -Inf) left.knot = tf.order

    # knots
    lk = left.knot#9
    tk = test.knot#20
    rk = right.knot#30

    # construct each part of v
    lv = seq(from=0,to=1,length=tk+1-lk)
    lv = lv[-length(lv)]
    rv = seq(from=1,to=0,length=rk+1-tk)
    ll = rep(0,lk-1)
    rr = rep(0,n-rk)
    cv = scale(c(lv,rv),center=TRUE,scale=FALSE)
    v = c(ll,cv,rr)
    
##    plot(v.direct,type='o')
#    plot(v.lrt,col='red',type='o')
#    lines(v/10, col = 'blue', type = 'o')
#    abline(v=final.model,col='lightgrey',lty=3)
    
    
    v = v/(max(abs(v))*10)
    v = v * (-sign.test.knot)
    
    return(v)
  }
  
  
  
  
  
  
  
  
  
  

# Function to aggregate successes and hits in order to calculate condit powers, for stopping time simulations
  getpowers.from.chunks = function(n, verdict.obj.name, ngrain=20, nsim = 100000, nchunks = 100, file.id = "bic-onejump-segmentsize-allbics"){

    # Helper functions
    getpow = function(verdicts,ii,loc){  return(sum(verdicts[ii,,loc],na.rm=T)/pmax(1,sum(!is.na(verdicts[ii,,loc]))))   }
    getpow.numer = function(verdicts){  return( sum(verdicts,na.rm=T))   }
    getpow.denom = function(verdicts){  return( pmax(1,sum(!is.na(verdicts))))   }

    # Powers
    powers = array(NA, c(ngrain,n))
    powers.numer = powers.denom = array(0,c(ngrain,n))
    powers.prox.numer = powers.prox.denom = powers.prox = rep(0,ngrain)
    
    loc = n/2
    proxlocs = (loc-log(n)):(loc+log(n))
    
    # extract numer and denom and sum over chunks
    for(chunk.i in 1:nchunks){
    cat('\r', chunk.i, "out of", nchunks)
    load(file=file.path(codedir, paste0("maxoutput/",file.id , n/2, "-chunk",chunk.i,".Rdata")))
#      verdict.obj.name = "verdicts.bic2"
      verdicts = eval(parse(text=verdict.obj.name))
      for(ii in 1:ngrain){
        # Exact powers
        for(loc in 1:n){
          powers.numer[ii,loc]  = powers.numer[ii,loc] + getpow.numer(verdicts[ii,,loc])
          powers.denom[ii,loc]  = powers.denom[ii,loc] + getpow.denom(verdicts[ii,,loc])
        }
        # Approximate powers at true break coordinate
        prox.verdicts = apply(verdicts[ii,,proxlocs], 1, 
                            function(this.sim.verdicts){
                             (if(!all(is.na(this.sim.verdicts))) any(this.sim.verdicts,na.rm=T) else NA)})
        powers.prox.numer[ii]  = powers.prox.numer[ii] + getpow.numer(prox.verdicts)
        powers.prox.denom[ii]  = powers.prox.denom[ii] + getpow.denom(prox.verdicts)
      }
      
      powers.prox.numer.bic = powers.prox.numer
      powers.prox.denom.bic = powers.prox.denom
    }
    cat('\n')

    # calculate condit power at each locaiton from extracted/summed numer and denom
    for(ii in 1:ngrain){
      powers[ii,]  = powers.numer[ii,]/powers.denom[ii,]
      powers.prox[ii] = powers.prox.numer[ii]/powers.prox.denom[ii]
    }
    rownames(powers) = names(powers.prox) = sigmalist
    colnames(powers) = 1:n

    return(list(powers=powers, powers.numer=powers.numer, powers.denom=powers.denom, powers.prox = powers.prox))
  }

# Function to calculate oracle powers
  getpowers.oracle = function(n,verdict.obj.name="verdicts.bic", ngrain=20, nsim = 100000, file.id = "bic-onejump-segmentsize-allbics"){
    # oracle power at true break coordinate
    powers.oracle = array(0, ngrain)
    for(chunk.i in 1:nchunks){
      cat('\r', chunk.i, "out of", nchunks)
      load(file=file.path(codedir, paste0("maxoutput/",file.id , n/2, "-chunk",chunk.i,".Rdata")))
      for(ii in 1:ngrain){
        powers.oracle[ii]  = powers.oracle[ii] + sum(verdicts.oracle[ii,],na.rm=T)/nsim
      }
    }
    names(powers.oracle) = sigmalist
    return(powers.oracle)
  }
  
# Function to extract stop times
# change stoptimes to conditional values (only the values for which loc=n/2 was tested)
  getstoptimes = function(n, ngrain=20, nsim = 100000, nchunks = 100, file.id = "bic-onejump-segmentsize-allbics"){
    loc = n/2
    stoptimes.bic.cond = stoptimes.bic2.cond = stoptimes.ebic.cond = stoptimes.sbic.cond = stoptimes.aic.cond = array(NA, c(ngrain,1))
    for(chunk.i in 1:nchunks){
      cat('\r', chunk.i, "out of", nchunks)
      load(file=file.path(codedir, paste0("maxoutput/",file.id , n/2, "-chunk",chunk.i,".Rdata")))
      # Collect conditional stoptimes in each chunk
      for(jj in 1:ngrain){ stoptimes.bic[jj,is.na(verdicts.bic[jj,,loc])] = NA }
      for(jj in 1:ngrain){ stoptimes.bic2[jj,is.na(verdicts.bic2[jj,,loc])] = NA }
      for(jj in 1:ngrain){ stoptimes.ebic[jj,is.na(verdicts.ebic[jj,,loc])] = NA }
      for(jj in 1:ngrain){ stoptimes.sbic[jj,is.na(verdicts.sbic[jj,,loc])] = NA }
      for(jj in 1:ngrain){ stoptimes.aic[jj,is.na(verdicts.aic[jj,,loc])] = NA }
      stoptimes.bic.cond = cbind(stoptimes.bic.cond,stoptimes.bic)
      stoptimes.bic2.cond = cbind(stoptimes.bic2.cond,stoptimes.bic2)
      stoptimes.ebic.cond = cbind(stoptimes.ebic.cond,stoptimes.ebic)
      stoptimes.sbic.cond = cbind(stoptimes.sbic.cond,stoptimes.sbic)
      stoptimes.aic.cond = cbind(stoptimes.aic.cond,stoptimes.aic)
    }
    return(list(stoptimes.bic = stoptimes.bic.cond, 
                stoptimes.bic2 = stoptimes.bic2.cond, 
                stoptimes.ebic = stoptimes.ebic.cond, 
                stoptimes.sbic = stoptimes.sbic.cond, 
                stoptimes.aic = stoptimes.aic.cond))
  }


## Function to make plot from 100 chunked simulations
makeplot = function(n, sigmalist, powers.list, stoptimes){
    #powers.list = powers.80    
    m = cbind(c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(1,1),c(2,2),c(2,2))
    layout(m)
    
    # Plot exact-location (n/2) condit powers
      loc = n/2
      cols = brewer.pal(8,'Dark2')
      xlim = c(0,max(sigmalist))
      plot(powers.list$powers.bic$powers[,loc] ~ sigmalist, type = 'l', lwd=2, ylim = c(0,1.5), xlim = xlim, axes = F, xlab="noise(sd)", ylab = "condit. powers", col = cols[1], cex.lab=2)
      title("Conditional Power at correct location, for one-jump", cex.main=2)
      axis(1, at = seq(from=0,to=3,by=0.5), label = seq(from=0,to=3,by=0.5))
      axis(2, at = seq(from=0,to=1,by=0.2), label = seq(from=0,to=1,by=0.2))
      lines(powers.list$powers.bic2$powers[,loc] ~ sigmalist, type = 'l', lwd=1, col = cols[2])
      lines(powers.list$powers.ebic$powers[,loc] ~ sigmalist, type = 'l', lwd=1, col = cols[3])
      lines(powers.list$powers.aic$powers[,loc] ~ sigmalist, type = 'l', col = cols[4], lwd=2)

    # Add proximal versions' powers
      lines(powers.list$powers.bic$powers.prox ~ sigmalist, type = 'l', lty=2, col = cols[1])
      lines(powers.list$powers.bic2$powers.prox ~ sigmalist, type = 'l', lty=2, col = cols[2])
      lines(powers.list$powers.ebic$powers.prox ~ sigmalist, type = 'l', lty=2, col = cols[3])
      lines(powers.list$powers.aic$powers.prox ~ sigmalist, type = 'l', col = 'red', lty=2)

    # Add fixed-stop powers
      lines(powers.list$powers.fixed1$powers[,loc] ~ sigmalist, type = 'l', col = 'green', lwd=2)
      lines(powers.list$powers.fixed2$powers[,loc] ~ sigmalist, type = 'l', col='darkgreen',lwd=2)

    # Plot oracle
      lines(powers.list$powers.oracle~sigmalist, col = cols[5])

    # Plot no-signal type II error
      abline(h = c(0.05,0.1,0.15,0.2), col = "grey", lty=2)
    
    # plot bic stoptimes
     addstoptimes = function(stoptimes, sigmalist, col, pch, adjust){
        par(new=TRUE)
        mn = apply(stoptimes,1,mean,na.rm=T) + 1 # adding one because of the way we define stoptime.
        sdev = apply(stoptimes,1,sd,na.rm=T)
        sigmalist_adjusted = sigmalist+rep(adjust,length(sigmalist))
        plot(mn ~ sigmalist_adjusted, pch=pch, ylim = c(0,7), xlim = xlim, axes=F, col = col, xlab = "", ylab = "")
        arrows(sigmalist_adjusted, mn-sdev, sigmalist_adjusted, mn+sdev, length=0.05, angle=90, code=3, col = col)

      }

      addstoptimes(stoptimes$stoptimes.bic,sigmalist, 'lightgrey', 19, 0)
      addstoptimes(stoptimes$stoptimes.bic2,sigmalist, 'lightgrey', "2", 0.05)
      addstoptimes(stoptimes$stoptimes.ebic,sigmalist, 'lightgrey', "e", 0.1)
      addstoptimes(stoptimes$stoptimes.aic,sigmalist, 'darkgrey', 19, -0.05)
      mtext("stoptimes", side = 4, padj = 6)
      axis(4, xlab="")
      
      # plot fixed stop times
      abline(h=c(1,2),lty = c(1,2), col = 'yellow')
      
    # Add Legend (on separate plot)
    plot(1,col='white',axes=F, xlab = "", ylab = "")
    legend("topleft", legend=c("BIC", "BIC-prox", "BIC-2", "BIC-2-prox", 
           "Extended BIC", "Extended BIC-prox", "AIC","AIC-prox", "Oracle"), 
           lwd = rep(1,11), lty = c(rep(1:2,5),1), col = c(rep(cols[1:4],each=2),cols[5]))
    # make legend
      legend("bottomleft", pch = c(19,19), col = c("lightgrey", "darkgrey"), 
             legend = c("avg bic stoptime+-1sdv","aic-stoptime"))
  }




