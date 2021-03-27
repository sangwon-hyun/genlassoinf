## Setup
library(genlassoinf)
library(dplyr)

## Generate data
sigma = 1
n = 100
maxsteps = 5
n = 300
sig = 2
beta = c(rep(0, n/3), rep(sig, n/3), rep(0, n/3))
set.seed(72) # use 13, 62 or 72
y = beta + rnorm(n, 0, sigma)

## Estimate fused lasso model
D = dual1d_Dmat(n)
f = dualpathSvd2(y, D, maxsteps,approx=T)

# Conduct naive & final model conditional inference on selected gaps
klater = 4
G = getGammat.naive(f, y, condition.step = klater)$G #k=klater,maxsteps=maxsteps,"dualpathSvd")
my.locs = f$action[1:klater]
u = rep(0, nrow(G))

## Produce p-values
pvals.spike = pvals.segment = pvals.naive = pvals.naive.ryan = c()
for(ii in 1:klater){
  d.spike = get_v_1dfusedlasso(f, y, k = ii, klater = klater, type = "spike")
  d.segment = get_v_1dfusedlasso(f, y, k = ii, klater = klater, type = "segment")
  pvals.spike[ii] = poly.pval(y, G, u, d.spike, sigma) %>% .$pv
  pvals.segment[ii] = poly.pval(y, G, u, d.segment, sigma) %>% .$pv
}

## Make plot
plot(y, pch=16, col='grey50', ylim = c(-3, 5))
lines(beta, col='red', lwd=3)
abline(v = f$cp[1:klater])
text(x=f$cp[1:klater], y=rep(4, klater) + (1:klater)/klater, 1:klater)

## Table of p-values
data.frame(loc = my.locs,
           spike = pvals.spike,
           segment = pvals.segment) %>% round(3)
