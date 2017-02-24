library(genlassoinf)

## Generate some data
consec = 2
D = makeDmat(60,type='trend.filtering',ord=0)
maxsteps = 10 
beta0 = rep(c(0,1),each=n/2)
set.seed(0)
y0 = beta0 + rnorm(n, 0, 1)
f0 = dualpathSvd2(y0, D, maxsteps, approx=T)
contrast = getdvec(f0,y0,test.step,stoptime.bic,type="segment")
ci.smuce    = contrast.smuce(y0,d)
