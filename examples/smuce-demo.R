library(stepR)
library(cghFLasso)

## First example: Generate SMUCE confidence bands, on a CGH example.
data(CGH)
y0 = CGH$GBM.y
sy = smuceR(y0, jumpint=T, confband=T)
cb = confband(sy)
plot(y0, col='grey50',pch=16)
lines(sy, lwd=2)
lines(cb, col='blue')


## Second example: simulated data
consec = 2
n = 60
D = makeDmat(n,type='tf',ord=0)
sigma = 1
lev1 = 0
lev2 = 1
maxsteps = 10 
beta0 = rep(c(lev1,lev2),each=n/2)
set.seed(1)
y0 = beta0 + rnorm(n, 0,sigma)


## Fit fused lasso path, and stop using BIC.
f0 = dualpathSvd2(y0,D,maxsteps,approx=T)
f0 = stop_path(f0, sigma=sigma, stoprule="bic")

    
## Identify test locations from path, stopping using BIC
locs = f0$pathobj$B[1:f0$stoptime]


## Collect SMUCE CIs and TG pvals, using the same segment contrasts
cis = array(NA,dim=c(n,2))
pvals = rep(NA,n)
for(test.step in 1:f0$stoptime){

    ## Form segment contrast for a particular location |loc|
    loc = locs[test.step]
    d = get_v_1dfusedlasso(f0,y0,test.step, f0$stoptime,type="segment")

    ## TG p-value
    pvals[loc] = poly.pval(y=y0, G=f0$Gobj.stoprule$G, u=f0$Gobj.stoprule$u, v=d, sigma=sigma)$pv

    ## SMUCE confidence interval
    cis[loc,1:2] = CI_smuce(y0,d)
}


## Optionally, create a step-sign plot /manually/; see step-sign-demo.R
signs = f0$ss[[f0$stoptime+1]]
final.model = f0$states[[f0$stoptime+1]]
s0 = step_sign_plot_inner(signs, final.model, n)# plot=TRUE
s0 = step_sign_plot(f0, stoptime = f0$stoptime, n, plot=TRUE)
