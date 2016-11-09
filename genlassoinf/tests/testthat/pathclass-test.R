context("Testing whether path class works well")

lev1 = 0 
lev2 = 2
sigma = 1
n = 10
beta0 = c(rep(lev1,n/2), rep(lev2,n/2))
set.seed(0)
y0 = beta0 + rnorm(n,0,sigma)
numsteps = 5
f0 <- dualpathSvd2(y0, dual1d_Dmat(length(y0)), maxsteps = numsteps,verbose=FALSE,approx=TRUE)

## Well, as of now, no real test here.
