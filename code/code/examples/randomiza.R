workingdir = "."
outputdir = "output/"

source('funs.R')
source('selectinf/selectiveInference/R/funs.inf.R')
source('dualPathSvd2.R')

sigma = .1
n = 100
maxsteps = 5
set.seed(72) # use 13, 62 or 72
nsim=10
klater=2
p.orig = p.new = matrix(NA,ncol=n, nrow = nsim)
for(isim in 1:nsim){
    lev1 = 0
    lev2 = 1
    beta0 = c(rep(lev1,n/2),rep(lev2,n/2))
    y0    = beta0 + rnorm(n,0,sigma)
    f0    = dualpathSvd2(y0,
                         D=dual1d_Dmat(length(y0)),
                         maxsteps,
                         approx=T)
   
    # Get orignal G
    G = getGammat.naive(f0,y0,condition.step=maxsteps)$G

    # Get my locations
    my.locs = f0$action[1:klater]
  
    ## Randomize and obtain stacked G
    newsigma=.01
    y0newlist = lapply(1:10, function(seed){set.seed(seed);rnorm(n,0,newsigma)})
    Gs = list()
    for(jj in 1:10){
        y0new = y0newlist[[jj]]
        f0new  = dualpathSvd2(y0new, D=dual1d_Dmat(length(y0new)), maxsteps,approx=T)
        Gs[[jj]] = getGammat.naive(f0new,y0new,condition.step=maxsteps)$G
    }
    newG = do.call(rbind,Gs)

    for(ii in 1:klater){
        d.segment = getdvec(f0,y0,k=ii,klater=klater,type="segment")
        
        ## Do inference
        p.orig[isim,my.locs[ii]] = pval.fl1d(y0,G,d.segment,sigma)
        p.new[isim,my.locs[ii]] = pval.fl1d(y0,newG,d.segment,sigma)
   }
}

head(p.orig)
p.orig[,48:52]
