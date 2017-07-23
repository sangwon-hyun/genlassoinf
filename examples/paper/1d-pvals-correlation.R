## MAKE sure you're working from [dropboxfolder]/code
source("settings.R")
source('selectinf/selectiveInference/R/funs.inf.R')
source('funs.R')
source('dualPathSvd2.R')
library(xtable)

# Make p-values
nsim=1000
maxsteps=5
lev1=0
lev2=0
sigma=1
n=60
D=dual1d_Dmat(n)
pseg = array(NA,dim=c(nsim,maxsteps-1))
for(isim in 1:nsim){
    print(isim)
    beta0 = rep(c(lev1,lev2),each=n/2)
    y0    = beta0 + rnorm(n, 0,sigma)
    f0    = dualpathSvd2(y=y0,D=D,maxsteps=maxsteps)
    for(istep in 1:(maxsteps-1)){
        G      = f0$Gammat[1:f0$nk[istep],]
        d      = getdvec(obj=f0, y=y0, k=istep, type="segment")
        my.pval = pval.fl1d(y=y0,G=G,d=d,sigma=sigma)
        pseg[isim,istep] = my.pval
    }
}

## make table
mytable = cor(pseg)
xtable(round(mytable,3))
