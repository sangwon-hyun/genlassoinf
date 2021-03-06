## Make sure you're working from [dropboxfolder]/code
source('funs.R')
source('testfuns.R')
source('pval.R')
source('dualPathSvd2.R')
source('selectinf/selectiveInference/R/funs.inf.R')
workingdir = "." 
outputdir = "output/"

n = 60
sigma = 1
lev1= 0
lev2= 2
lev2list=c(0,1,2)
numsteps=1

set.seed(0)
sigma = 2
lev1= 0
pt.cex=.5  
beta0 = c(rep(0,20),rep(5,10),rep(10,25)) + 10
set.seed(33)
y0  = beta0 + rnorm(n=length(beta0),sd=sigma)

maxsteps = 5
k = 2
f0  = dualpathSvd2(y0, D=dual1d_Dmat(length(y0)), maxsteps,approx=T)  
G = getGammat.naive(f0,y0,condition.step=maxsteps)$G
pvals.spike = pvals.segment = c()
my.locs = f0$action[1:k]
my.ord = order(my.locs)
my.locs = sort(my.locs)

for (i in 1:k){
  d.spike = getdvec(f0,y0,k=i,klater=k,type="spike")
  d.segment = getdvec(f0,y0,k=i,klater=k,type="segment")
  pvals.spike[i] = pval.fl1d(y0,G,d.spike,sigma)
  pvals.segment[i] = pval.fl1d(y0,G,d.segment,sigma)
}

my.pvals.spike = pvals.spike[my.ord]
my.pvals.segment = pvals.segment[my.ord]

xlab = "Location"
w = 5; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-8,30)
mar = c(4.5,4.5,0.5,0.5)
let = c("A","B")
col.hline='lightgrey'

# Plot, data and segment contrast
pdf(file=file.path(outputdir, "1dfl-segment.pdf"),width=w,height=h)
par(mar=mar)
plot(y0, ylim = ylim,axes=F, xlab = xlab, ylab = "", pch=pch, col=pcol);
lines(f0$beta[,k+1],col="blue",lwd=lwd)
axis(1); axis(2)
abline(h=0,col=col.hline)
for (i in 1:k) {
  Ji = my.locs[i]
  JL = ifelse(i==1, 0, my.locs[i-1])
  JR = ifelse(i==k, n, my.locs[i+1])
  v.segment = rep(0,n)
  v.segment[(JL+1):Ji] = -1/(Ji-JL)
  v.segment[(Ji+1):JR] = 1/(JR-Ji)
  if (i==2) points(v.segment*50, pch=17, col=3) # only draw the 2nd one
  text(let[i], x=Ji+2, y=30)
  text(sprintf("p-value = %0.3f",my.pvals.segment[i]),
       pos=ifelse(i==1,2,4), x=Ji, y=8-(k-i)*2)
}
abline(v=my.locs, lty=2)
legend("topleft", pch=c(pch,NA,NA,17), lty=c(NA,1,2,NA), lwd=c(NA,lwd,1,NA),
       col = c(pcol,"blue","black",3),legend=c("Data","Estimate",
                                        "Changepoint","Segment contrast"))
graphics.off()

# Plot, data and spike contrast
pdf(file=file.path(outputdir, "1dfl-spike.pdf"),width=w,height=h)
par(mar=mar)
plot(y0, ylim = ylim,axes=F, xlab = xlab, ylab = "", pch=pch, col=pcol);
lines(f0$beta[,k+1],col="blue",lwd=lwd)
axis(1); axis(2)
abline(h=0,col=col.hline)
for (i in 1:k) {
  Ji = my.locs[i]
  v.spike = rep(0,n)
  v.spike[Ji] = -1
  v.spike[Ji+1] = 1
  if (i==2) points(v.spike*5, pch=17, col=2) # only draw the 2nd one
  text(let[i], x=Ji+2, y=30)
  text(sprintf("p-value = %0.3f",my.pvals.spike[i]),
       pos=ifelse(i==1,2,4), x=Ji, y=8-(k-i)*2)
}
abline(v=my.locs, lty=2)
legend("topleft", pch=c(pch,NA,NA,17), lty=c(NA,1,2,NA), lwd=c(NA,lwd,1,NA),
       col = c(pcol,"blue","black",2),legend=c("Data","Estimate",
                                        "Changepoint","Spike contrast"))
graphics.off()

####################################################
## ### 1-jump Example ######

##     # generate stuff
##     nsim = 10000
##     n = 60
##     sigma = 1
##     lev1= 0
##     lev2= 2
##     lev2list=c(0,1,2)
##     numsteps=1
##     spike   = introexample(testtype = "spike", nsim=nsim,sigma=sigma,lev1=lev1,lev2=lev2,lev2list=lev2list,numsteps=numsteps,verbose=T)
##     segment = introexample(testtype = "segment", nsim=nsim,sigma=sigma,lev1=lev1,lev2=lev2,lev2list=lev2list,numsteps=numsteps,verbose=T)
##   #save(file = file.path(outputdir,"onejump-example.Rdata"), list = c("spike","segment","proportions0", "proportions1"))
##   load(file.path(outputdir,"onejump-example.Rdata"))

##   pdf(file.path(outputdir,"spiketest-example.pdf"),width=10,height=5)
##     par(mfrow = c(1,2))
##     set.seed(0)
##     sigma = 2
##     lev1= 0
##     pt.cex=.5
##     ylim = c(-10,25)
##     #beta0 = c(rep(5,10),rep(3,10),rep(15,10),rep(0,10),rep(5,10),rep(0,10))
##     beta0 = c(rep(0,20),rep(5,10),rep(6,15),rep(7,15))
##     set.seed(2)
##     y0    = beta0 + rnorm(n=length(beta0),sd=sigma)
##     plot(y0, axes=F, cex=pt.cex, ylim = ylim, xlab = 'coordinate', ylab = 'y')
##     lines(beta0, col='red',lwd=1.5)
##     abline(v = 20, lty=2, col='lightgrey')
##     text(x=21,y=-8, labels=expression(i[k]))
##     axis(1); axis(2);
##     title(main = expression(paste("Data (y =",  theta + epsilon, ") and Test Location (", i[k],")")))
##     legend("topright", pch = c(1,NA), pt.cex=pt.cex, lty = c(NA,2), col = c('black', 'lightgrey'), legend = c("data", "test location"))
    
##     plot(x=c(20,21),y=c(-1,1), col = 'blue', pch = 17, xlim = c(0,60), ylim = c(-2,4), axes=F, xlab = 'coordinate', ylab = 'y')
##     abline(v = 20, lty=2, col='lightgrey')
##     points(x = c(1:19, NA,NA,22:60),y = c(rep(0,19),NA,NA,rep(0,39)), pch = 17, col = 'blue')
##     points(x = c(1:20,21:60),y = c(rep(-0.3,20),rep(0.3,40)), col = 'green', pch = 16)
##     text(x=19,y=4, labels=expression(i[k]))
##     axis(1); axis(2);
##     title(main = expression(paste("Spike and Segment Test Contrasts (", d[i[k]],")")))
##     legend("topright", pch = c(17,16, NA), lty = c(NA,NA,2), col = c('blue',  'green','lightgrey'), legend = c("spike contrasts","segment contrasts", "test location"))

##     dev.off()
     
