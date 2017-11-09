## Synopsis: Make an example plot of a one-sided confidence interval, in a 1d, one-jump example.
library(binSegInf)
outputdir = "../output"

## Run things once.
n=60
set.seed(53)
sigma=1
lev=3
mn = c(rep(0,n/2), rep(lev,n/2))
y = mn + rnorm(n,0,sigma)
g = dualpathSvd2(y, maxsteps=1, D=makeDmat(n,ord=0,type='tf'))
p = polyhedra_from_genlasso(g)
contrast = make_all_segment_contrasts(g)[[1]]
spike.contrast = rep(0,n)
spike.contrast[g$cp] = -g$cp.sign
spike.contrast[g$cp+1] = +g$cp.sign
alpha=.05
myci = confidence_interval(y,p, contrast,sigma=1,alpha=alpha,alternative="one.sided", fac=30, gridsize=10000)
myci.spike = confidence_interval(y,p, spike.contrast,sigma=1,alpha=alpha,alternative="one.sided", fac=30, gridsize=10000)

## Calculate a few useful things.
left.sample.mean = mean(y[1:g$cp[1]])
right.sample.mean = mean(y[(g$cp[1]+1):n])
middle = mean(c(left.sample.mean, right.sample.mean))
rng = c(middle - myci[1]/2, middle + myci[1]/2)
rng.spike = c(middle - myci.spike[1]/2, middle + myci.spike[1]/2)


## Plot settings
pcol.sample.mean = "coral1"
lcol.middle = "grey70"
lty.sample.mean = "solid"
lty.middle = lty.ci = "dashed"
pcol.data = "grey50"
cols = RColorBrewer::brewer.pal(8,"Set2")
lcol.ci = cols[3]
lcol.ci.spike = cols[4]
pch.data = 16
xlab = "Location"
lwd.ci = 2.5
lwd.means = 2.5
lwd.arrows = 2.5

## Plot it:
filename = "ci-example.pdf"
## pdf(file=file.path(outputdir, filename), width=5, height=5)
par(mar=c(4.1,2.1,0.6,0.5))
plot(y, ylim=range(y)+c(-0.5,0), xlim = c(0,65), col=pcol.data, pch=pch.data, axes=FALSE, ylab="", xlab = xlab)
axis(1); axis(2)
lines(x=1:(g$cp[1]), y=rep(left.sample.mean, g$cp[1]), col=pcol.sample.mean, lty = lty.sample.mean, lwd=lwd.means)
lines(x=(g$cp[1]+1):n, y=rep(right.sample.mean, n-g$cp[1]), col=pcol.sample.mean, lty = lty.sample.mean, lwd=lwd.means)
abline(h=middle,col=lcol.middle, lty=lty.middle)
abline(h=rng, col = lcol.ci, lwd=lwd.ci, lty=lty.ci)
abline(h=rng.spike, col = lcol.ci.spike, lwd=lwd.ci, lty=lty.ci)
arrows(x0=g$cp[1]+0.5,
       x1=g$cp[1]+0.5,
       y0=rng[1],
       y1=rng[2],
       code=3,
       length=0.1,
       lwd=lwd.arrows)

arrows(x0=g$cp[1] - 1,
       x1=g$cp[1] - 1,
       y0=rng.spike[1],
       y1=rng.spike[2],
       code=3,
       length=0.1,
       lwd=lwd.arrows)
text(label=bquote(eta[0.05]^spike),x=g$cp[1]-4, y=middle-.2, cex=1.1)
text(label=bquote(eta[0.05]^segment),x=g$cp[1]+6, y=middle-0.9, cex=1.2)

legend("bottomright",
       legend=c("Sample means","Midpoint", "Data", "Spike CI lower bound", "Segment CI lower bound"),
       lty=c(lty.sample.mean, lty.middle, NA, lty.ci, lty.ci),
       col=c(pcol.sample.mean, lcol.middle, pcol.data, lcol.ci.spike, lcol.ci),
       lwd=c(lwd.means,1, NA, lwd.ci, lwd.ci),
       pch = c(NA,NA, pch.data,NA,NA))
graphics.off()
