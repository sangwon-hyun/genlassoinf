# Generate data
sigma = .1
n = 100
maxsteps = 5
lev1 = 0
lev2 = 2
consec = 2
D = makeDmat(n,type='tf',ord=0)

## Generate data + path
beta0 = rep(c(lev1,lev2),each=n/2)
set.seed(0)
y0    = beta0 + rnorm(n, 0,sigma)
f0    = dualpathSvd2(y0, D=D, maxsteps, approx=T)


## Method 1: Get naive poyhedron (fixed stop time of 1)
states = get.states(f0$action)
stoptime = 1
Gobj.naive = getGammat.naive(obj = f0, y = y0, condition.step = 1)
G = Gobj.naive$G
u = Gobj.naive$u


## Method 2: First way to /manually/ get stop-time-incorporated polyhedron.
mm = get.modelinfo(obj=f0, consec=2, sigma=sigma)
bic = mm$ic
stoptime = which.rise(bic,consec) - 1
stoptime = pmin(stoptime, n-consec-1)
Gobj.stoprule = getGammat.with.stoprule(obj = f0, y = y0,
                                        condition.step = stoptime+consec, type ='tf',
                                        stoprule = "bic", sigma = sigma,
                                        consec = consec, maxsteps = maxsteps, D = D)
G = Gobj.stoprule$G
u = Gobj.stoprule$u

## Method 3: Better (more user-friendly) way to stop using BIC, add the information to f0,
## then extract the polyhedron
f0 = stop_path(f0, sigma=sigma, stoprule="bic")
G = f0$Gobj.stoprule$G
u = f0$Gobj.stoprule$u
stoptime = f0$stoptime



## AFTER running ONE OF Method 1-3: Form contrast and segment test p-value
locs = f0$stoppedmodel
pvals = rep(NA,length(locs))
vs = list()
for(ii in 1:length(locs)){
    this.sign = f0$pathobj$s[which(f0$pathobj$B == locs[ii])]
    my.v.lrt = make.v.tf.fp(test.knot = locs[ii],
                            adj.knot  = locs,
                            test.knot.sign = this.sign,
                            D=D)
    pval = poly.pval(y=y0, G=G, u=u, v=my.v.lrt, sigma=sigma)$pv
    cat("After fixed size model was selected, the segment test pvalue at location",
        locs[ii], "is", pval,fill=TRUE)
    pvals[ii] = pval
    vs[[ii]] = my.v.lrt
}

## Plot the p-values on the data.
xlab = "Location"
w = 5; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-3,5)
mar = c(4.5,4.5,0.5,0.5)
xlim = c(0,110)

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
pcols.delta =   pcols.oneoff = RColorBrewer::brewer.pal(n=3,name="Set2")
pch.spike = 15
pch.segment = 17
cex.contrast = 1.2

par(mar=c(4.1,3.1,3.1,1.1))
plot(y0, ylim = ylim,axes=F, xlim=xlim, xlab = xlab, ylab = "", pch=pch, col=pcol);
axis(1);axis(2)
abline(v=f0$stoppedmodel, col='lightgrey', lwd=2, lty=3)
text(x=f0$stoppedmodel,y=4,label=paste0("segment test p-value = ", pval))

## CI
myci = ci(y0, G, u, vs[[1]], sigma = sigma, alpha = 0.95, alternative="two.sided")
myci


mabline(h = mean(v.segment,na.rm=T), col = 'lightgrey')
legend("topleft", pch=c(pch.dat,NA,pch.spike,pch.segment),
       lty=c(NA,1,NA,NA), lwd=c(NA,2,NA,NA),
       col = c(pcol.dat, lcol.sig, pcol.spike,pcol.segment),
       pt.cex = c(cex.contrast, NA, cex.contrast, cex.contrast),
       legend=c("Data", "Mean","Spike contrast", "Segment contrast"))
title(main=expression("Data example"))
graphics.off()



## Also produce confidence intervals and plot them.

