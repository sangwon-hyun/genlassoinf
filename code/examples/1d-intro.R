## Make sure you're working from [dropboxfolder]/code
# Ryan
workingdir = "."
outputdir = "output/"

#source('helpers.R') # generic helpers
source('funs.R')
source('testfuns.R')
source('selectinf/selectiveInference/R/funs.inf.R')
source('dualPathSvd2.R')
library(xtable)

# Helper Function
pval.naive = function(y0,loc,sigma){
  gap = y0[loc+1] - y0[loc]
  positive.gap = sign(gap) * gap
  #threshold = qnorm(0.975,mean=0,sd=sigma)
  pval = 2*(1-pnorm(positive.gap, mean=0, sd = sqrt(2*sigma^2)))
  return(pval)
}

pval.naive.ryan = function(y0, loc, loc.left, loc.right, sigma) {
  gap = abs(mean(y0[loc.left:(loc-1)]) - mean(y0[loc:(loc.right-1)]))
  sd = sqrt(sigma^2*(1/(loc-loc.left)+1/(loc.right-loc)))
  return(2*(1-pnorm(gap,mean=0,sd=sd)))
}

# Visualize
sigma = .1
n = 100
maxsteps = 5
#  for(myseed in 10:100){
#    readline()
#    set.seed(myseed)
    
set.seed(72) # use 13, 62 or 72
beta0 = onejump.y(returnbeta=T, lev1=0, lev2=1, sigma, n)
y0    = onejump.y(returnbeta=F, lev1=0, lev2=1, sigma, n)
f0    = dualpathSvd2(y0, D=dual1d_Dmat(length(y0)), maxsteps,approx=T)
   
# Conduct naive & final model conditional inference on selected gaps
klater = 2
G = getGammat.naive(f0,y0,condition.step=klater)$G#k=klater,maxsteps=maxsteps,"dualpathSvd")
pvals.spike = pvals.segment = pvals.naive = pvals.naive.ryan = c()
my.locs = f0$action[1:klater]
my.ord = order(my.locs)
my.locs = sort(my.locs)

for(ii in 1:klater){
  d.spike = getdvec(f0,y0,k=ii,klater=klater,type="spike")
  d.segment = getdvec(f0,y0,k=ii,klater=klater,type="segment")
                                        #stopifnot(d[f0$action[ii]]!=0)
  pvals.spike[ii] = pval.fl1d(y0,G,d.spike,sigma)
  pvals.segment[ii] = pval.fl1d(y0,G,d.segment,sigma)
  pvals.naive[ii] = pval.naive(y0,loc=f0$action[ii],sigma=sigma)
  pvals.naive.ryan[my.ord[ii]] = pval.naive.ryan(y0,loc=my.locs[ii],
                    loc.left=ifelse(ii==1,1,my.locs[ii-1]),
                    loc.right=ifelse(ii==klater,n+1,my.locs[ii+1]),
                    sigma=sigma)
}

cols.naive = (pvals.naive < 0.05/klater)+1
cols.spike = (pvals.spike < 0.05/klater)+1
cols.segment = (pvals.segment < 0.05/klater)+1
  
xlab = "Location"
w = 5; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-0.5,1.5)
mar = c(4.5,4.5,0.5,0.5)

# Plot, data and mean
pdf(file=file.path(outputdir, "intro1.pdf"),width=w,height=h)
par(mar=mar)
plot(y0, ylim = ylim,axes=F, xlab = xlab, ylab = "", pch=pch, col=pcol);
lines(beta0,col='red',lwd=lwd)
axis(1); axis(2)
legend("bottomright", pch=c(pch,NA), lty=c(NA,1), lwd=c(NA,lwd),
       col = c(pcol,"red"),legend=c("Data","Mean"))
graphics.off()

# Plot, data and estimate
pdf(file=file.path(outputdir, "intro2.pdf"),width=w,height=h)
par(mar=mar)
plot(y0, ylim = ylim,axes=F, xlab = xlab, ylab = "", pch=pch, col=pcol);
lines(f0$beta[,klater+1],col='blue',lwd=lwd)
axis(1); axis(2)
text(c("A","B"), x = my.locs+5, y = rep(1.5,1.5))
abline(v = my.locs, col = "black", lty=2)
legend("bottomright", pch=c(pch,NA,NA), lty=c(NA,1,2),lwd=c(NA,lwd,1),
       col = c(pcol,"blue","black"),legend=c("Data","Estimate","Changepoint"))
graphics.off()

# Generate table of p-values
un.ord = order(my.ord)
mytable = cbind(my.locs,pvals.naive.ryan[un.ord],pvals.segment[un.ord])
rownames(mytable) = c("A","B")
colnames(mytable) = c("Location","Naive p-values","TG p-values")
xtable(mytable, digits=c(1,0,3,3), align="|r|r|r|r|")





# Collect the true jump-points's naive p-value distribution, and null jump-point's naive p-value distribution'
  pvals.null = pvals.nonnull = c()
  for(isim in 1:1000){
    beta0 = onejump.y(returnbeta=T, lev1=0, lev2=1, sigma, n)
    y0    = onejump.y(returnbeta=F, lev1=0, lev2=1, sigma, n)

    pvals.null[isim] = pval.naive(y0,48,sigma)
    pvals.nonnull[isim] = pval.naive(y0,50,sigma)
  }
  par(mfrow=c(1,2))
  hist(pvals.null,xlim = c(0,1))
  hist(pvals.nonnull, xlim=c(0,1))  
