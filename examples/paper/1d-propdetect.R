source('funs.R')
library(genlasso)

## Simulations settings 
lev1=0
lev2=0
maxsteps=1
n  = 60
sigma=1
mysim = function(lev1, lev2){
    beta0 = rep(c(lev1,lev2),each=n/2)
    y0    = beta0 + rnorm(n, 0,sigma)
    D = makeDmat(n, ord=0)
    f0    = genlasso(y0,D=D,maxsteps=maxsteps,approx=T)
    return(f0$pathobj$B)
}

## Get proportions of detections
props = rep(NA,3)
nsim = 1000
lev2s = c(0,1,2)
for(lev2 in lev2s){
    splits = replicate(nsim, mysim(0,lev2))
    hist(splits)
    props[which(lev2==lev2s)] = sum(splits==30)/nsim
}
names(props) = paste("gap=",lev2s)
print(props)
