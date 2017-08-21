## MAKE sure you're working from [dropboxfolder]/code
## source("settings.R")
## source('selectinf/selectiveInference/R/funs.inf.R')
## source('funs.R')
## source('dualPathSvd2.R')
## library(genlasso)
## library(Matrix)
## library(RColorBrewer)
load(file="../main/paper/data/stockData.RData")
outputdir = "../main/paper/data/"

# Format data and response
center = function(vec){
  as.numeric(scale(vec, center=TRUE,scale=FALSE))
}
log.serial.transform = function(vec){
  ratios =  vec[2:length(vec)] / vec[1:(length(vec)-1)]
  logratios = log(ratios)
  return(logratios)
}
X.uncentered = tricky_prices[,c(1,2,3)]#c(1,2,3,4,5)
## matplot(X.uncentered)
## matplot(tricky_prices[,30],type='o',pch=sapply(1:30,toString))
X = apply(X.uncentered,2,log.serial.transform)
X = apply(X,2,center)
TT = nrow(X)
J = ncol(X)

## This setting with any noise detects a spurious cut
beta0 = c(rep(c(-10,10,-10),each=TT/3),-10,-10,
          ## rep(c(-10,10),each=TT/2),10,  # (-5 and 5)
          ## rep(c(10,10,10),each=TT))   #
          rep(10,each=TT))/10

load(file=file.path(outputdir, paste0("stockData-sim.Rdata")))

X.augmented = do.call(cbind, lapply(1:ncol(X), function(icol) diag(X[,icol])))
mn = X.augmented %*% matrix(beta0,nrow=length(beta0))

betalist.ind = lapply(1:5, function(ii) (TT*(ii-1)+1):(TT*ii))
membership   = sapply(final.model.decluttered,
                      function(thiscut) which(unlist(lapply(betalist.ind,
                                                            function(thisvar.ind) any(thiscut %in% thisvar.ind)))))
membership.before.declutter = sapply(final.model.orig,
                                     function(thiscut) which(unlist(lapply(betalist.ind,
                                                                           function(thisvar.ind) any(thiscut %in% thisvar.ind)))))


#############
## Plot  ####
#############
ylim = c(0,20)
xlim = c(0,TT+10)
x = c(1:TT)
pch=16
w=6;h=3;
xlab = "Location"
ylab = ""
## lty.knot=3
lty.knot=2
## lwd.knot=1.5
lcol.knot="blue"#"grey40"
## lcol.knot.before.declutter="grey85" ## It was pointed out that this was too light
## lcol.knot.before.declutter="grey50"
cols = lcols.beta = RColorBrewer::brewer.pal(3,"Set2")

lwd.beta = 2
lty.beta = 1
betalist = list(beta0[1:TT],
                beta0[(TT+1):(2*TT)],
                beta0[(2*TT+1):(3*TT)])
xx = (final.model.decluttered - (membership-1)*TT)
xx.before.declutter = (final.model.orig - (membership.before.declutter-1)*TT)
mar = c(2,3.5,1,1)
mar3 = c(3.5,3.5,1,1)
my.est = MASS::ginv(X.tilde) %*% f0$beta[,stop.time+1]
my.est.list = list(my.est[1:TT],
                   my.est[(TT+1):(2*TT)],
                   my.est[(2*TT+1):(3*TT)])
lcol.est = 'black'
lwd.est = 2
lty.est = 1
xlim=c(0,TT+50)
jjth = c("1st","2nd","3rd")

for(jj in 1:3){
    if(jj==3){ h=2.6 }
    if(jj==2){ h=2.6 }
    if(jj==1){ h=2 }
    ylab = paste(jjth[jj], "coefficient")
    pdf(file.path(outputdir,paste0("new-stockdata-",jj,".pdf")),width=w,height=h)

    if(jj==3){ par(mar=mar3)} else {par(mar=mar)}

    ## collect some things
    rng = range(betalist[[jj]])
    this.xx = xx[membership==jj]
    this.xx.before.declutter = xx.before.declutter[membership.before.declutter==jj]
    this.beta = betalist[[jj]]
    this.pseg = signif(pseg[membership==jj],3)

    ## draw signal and detected knots
    ylim = rng + (rng[2]-rng[1])/3*c(-1,1)
    if(jj==1){ ylim[2] = ylim[2] + 0.8; ylim[1] = ylim[1]+0.2}
    ## xlim = c(1, (TT+10))
    plot(NA,ylim=ylim,xlim=xlim,axes=F,xlab="",ylab="")
    abline(v = this.xx.before.declutter, lty=lty.knot, lwd=lwd.knot, col = lcol.knot.before.declutter)
    abline(v = this.xx,lty=lty.knot, lwd=lwd.knot,col=lcol.knot)

    lines(x=(1:TT)-jj+1, y=this.beta, col = lcols.beta[jj], lwd = lwd.beta)
    if(jj==3){ axis(1);title(xlab=xlab,line=1.9) } else { axis(1)}
    axis(2);{title(ylab=ylab, line = 2.1)}

    ## Label locations
    myletters = toupper(letters[1:length(this.xx)] )
    this.yy = ylim[2]-(ylim[2]-ylim[1])/10 + (-sum(membership==jj)+1):(0)/10
    if(jj==3) this.yy = this.yy+0.2
    text(x=this.xx+5, y=this.yy, label = myletters)

    position = (if(jj==1) "topright" else "bottomright")

    ## Add the estimate lines
    lines(x=(1:TT)-jj+1, y=my.est.list[[jj]], lwd = lwd.est, col = lcol.est)

    inset=0

    ## Show p-values in legend
    if(!all(is.na(this.pseg))){
        if(jj==2){
        legend(position,
               lty = c(rep(NA,length(myletters)),rep(lty.knot,2),lty.est,lty.beta),
               lwd = c(rep(NA,length(myletters)),rep(lwd.knot,2),lwd.est,lwd.beta),
               col = c(rep("black",length(myletters)),lcol.knot, lcol.knot.before.declutter,lcol.est,lcols.beta[jj]),
               pch = c(myletters, rep(NA,4)),
               legend = c(paste("p-value =",this.pseg),
                 "Decluttered changepoints","Original changepoints", "Estimate","Truth"),
                             ## paste(myletters,
                             ##  rep(":", length(this.pseg)),
               ## pch=rep(NA,length(this.pseg)),
               ## title="p-values by location",
               bg="white",
               x.intersp=1,
               bty="o",
               box.lwd = 1)
        ## text.width=c(0.5,0.5,0.5,0.5))
       } else {
        legend(position,
               lty = c(rep(NA,length(myletters)),lty.est,lty.beta),
               lwd = c(rep(NA,length(myletters)),lwd.est,lwd.beta),
               col = c(rep("black",length(myletters)),lcol.est,lcols.beta[jj]),
               pch = c(myletters,rep(NA,2)),
               legend = c(paste("p-value =",this.pseg),"Estimate", "Truth"),
                             ## paste(myletters,
                             ##  rep(":", length(this.pseg)),
               ## pch=rep(NA,length(this.pseg)),
               ## title="p-values by location",
               bg="white",
               x.intersp=1,
               bty="o",
               box.lwd = 1, inset=inset)
               ## text.width=c(0.5,0.5,0.5,0.5))
   }
    }
   graphics.off()
}


################################
## Plot just the five stocks ###
################################
ltys = rep(1,5)
xlab="Trading day (2015)"
ylab="Stock price"
xlim=c(0,270)
w=5;h=7.5
mar=c(4,1,1,1)

pdf(file=file.path(outputdir,"stockdata-raw.pdf"),width=w,height=h)
matplot(tricky_prices[,1:3],col=cols,type='l',lty=ltys, ylim=c(60,190),axes=F,ylab=ylab,xlab=xlab)
text(x=rep(260,3)-20,
     y=as.numeric(tricky_prices[TT-1,1:3])+c(10,5,15),
     labels=paste(rep("stock",3),1:3))
## lines(scale(exp(y[2:TT]/y[1:(TT-1)]))*10+100,type='l')
## klines(y+100,col=cols,type='l')
axis(1)
axis(2)
graphics.off()
