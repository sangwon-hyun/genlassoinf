library(stepR)
library(cghFLasso)

## Generate some data
consec = 2
n = 60
D = makeDmat(n,type='tf',ord=0)
sigma = 1
lev1 = 0
lev2 = 3
beta0 = rep(c(lev1,lev2,lev1),each=n/3)
set.seed(35)
y0 = beta0 + rnorm(n, 0,sigma)

## Fit fused lasso path, and stop using BIC.
f0 = dualpathSvd2(y0,D,10,approx=T)
f0 = stop_path(f0, sigma=sigma, stoprule="bic")

## Produce step-sign plot manually
signs = f0$ss[[f0$stoptime+1]]
final.model = f0$states[[f0$stoptime+1]]
s0 = step_sign_plot_inner(signs, final.model, n)# plot=TRUE

## Or do it with the path object |f0| alone!
s0 = step_sign_plot(f0, f0$stoptime, n, plot=TRUE)

## From this plot, the user (YOU!!) should make decision of what segment
## contrast to test. You can opt to use an automatic post-processing rule for
## breakpoints, such as centroid clustering.

## Based on the plot from |s0|, let's say you decided that out of c(18,20,40),
## 18 and 20 are too closeby, so of the two, you want to test only 20, so that
## c(20,40) are to be tested, with mean shift directions (+1,-1). Manually form
## the two segment contrasts and conduct TG tests

d = make_contrast(20,c(1,40),+1,60)
poly.pval(y=y0, G=f0$Gobj.stoprule$G, u=f0$Gobj.stoprule$u, v=d, sigma=sigma)$pv

d = make_contrast(40,c(20,60),-1,60)
poly.pval(y=y0, G=f0$Gobj.stoprule$G, u=f0$Gobj.stoprule$u, v=d, sigma=sigma)$pv
