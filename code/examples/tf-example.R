  # Make sure you're working from [dropboxfolder]/code
  source('funs.R')
  source('testfuns.R')
  source('dualPathSvd2.R')
  source('selectinf/selectiveInference/R/funs.inf.R')
  library(genlasso)
  library(RColorBrewer)
  verbose = F

  workingdir = "."
  outputdir = "output/"

#########################################
## Generate p-values for simulation  ####
#########################################

# Picture of data and signal
  beta1 = seq(from=1,to=20,by=1)
  beta2 = -0.5*seq(from=21,to=40,by=1) + 30
  beta3 = seq(from=41,to=60,by=1) - 30
  beta0 = c(beta1,beta2,beta3)

  bic.fac=1
  myseed=3#myseed = 42
  set.seed(myseed)
  sigma=5
  y0 = beta0 + rnorm(length(beta0),0,sigma)
  n = length(y0)

  maxsteps = 3
  stop.time = 2
  tf.order = 1
  D = makeDmat(n,type="trend.filtering",order=1)
  f0 = dualpathSvd2(y0, D, maxsteps = maxsteps)
  states = get.states(f0$action)
  final.model = states[[stop.time+1]]
  model.order = order(final.model)
  final.model = final.model[model.order]
  signs = f0$pathobj$s[model.order]

  #  ii = which(final.model==40); 
  pval.direct = pval.lrt = rep(NA,2)
  fits = v.lrt = v.direct =list()
  for(ii in 1:2){

    f1 = dualpathSvd2(y0,D,stop.time)
    this.sign = f1$pathobj$s[which(f1$pathobj$B == final.model[ii])]
    my.v.lrt = make.v.tf.fp(test.knot = final.model[ii],
                        adj.knot  = final.model,
                        test.knot.sign = this.sign, D=D)
    # get p-value
    Gobj.naive = getGammat.naive(obj=f0,y=y0,
                                 condition.step = stop.time)
    G = Gobj.naive$G
    u = rep(0,nrow(G))
    v.lrt[[ii]] = my.v.lrt
    mypval = pval.fl1d(y0, G, dik = v.lrt[[ii]],    sigma, u = u)
    print(mypval)
    pval.lrt[ii]     = signif(mypval,3)
    f1 = dualpathSvd2(y0,D,maxsteps=3)
    fits[[ii]] = f1$beta[,3]
  }
  
  
    
  ###################################
  # Plot trend filtering example ####
  ###################################

  # plot parameters
  xlab = "Location"
  ylab1 = ""#"Data and estimate"
  ylab2 = ""#"contrast(=v)"
  w = 5; h = 5
  pch.dat = 16; 
  pch.lrt = 17
  pcol.dat = "gray50"
  pcol.lrt = brewer.pal(n=3,name="Set2")[1]
  lcol.fit = "blue"
  lcol.adj.knot = lcol.hline = "lightgrey"
  lcol.test.knot = "blue"
  lty.fit = 1 ; lwd.fit = 2
  lty.knot = 2 ; lwd.knot = 1

  w = 5; h = 5
  ylim = c(-8,30)
  mar = c(4.5,4.5,0.5,0.5)
  rlab = seq(from=-1,to=1,by=0.1)
  rlab.pos = 50*(seq(from=-1,to=1,by=0.1))
  rlab[14:20] = rlab.pos[14:20]= NA

  for(ii in 1:2){
    pdf(file=file.path(outputdir, paste0("tf-examples-",ii,".pdf")),width=w,height=h)  
    par(mar=mar)
    plot(y0, ylim = c(-13,45), axes=F, pch = pch.dat, col=pcol.dat, xlab = xlab,ylab=ylab1); 
    lines(fits[[ii]], col = 'blue', lwd = lwd.fit)
    abline(v=final.model+1, lty=lty.knot, lwd=lwd.knot)
    points(v.lrt[[ii]]*30,col=pcol.lrt,pch=pch.lrt);
    abline(h=0, col = lcol.hline)
    text(x=final.model[ii]-1, y = -11, label = paste("p-value =", round(pval.lrt[ii],3)))
    text(c("A","B"), x=final.model+3, y=28)
    axis(1);axis(2);
    legend("topleft", legend=c("Data","Estimate","Knot","Segment contrast"),
            pch=c(pch.dat,NA,NA,pch.lrt),
            lty=c(NA,lty.fit,lty.knot,NA),
            lwd=c(NA,lwd.fit,lwd.knot,NA),
            col=c(pcol.dat,lcol.fit,1,pcol.lrt), bg="white")
    graphics.off()
  }
