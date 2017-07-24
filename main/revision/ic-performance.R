## Generate data
source("../main/revision/ic-performance-helper.R")
sigma = 1
n = 60
consec = 2

onesim <- function(lev,n,fun=onejump,isim=NULL){
    if(!is.null(isim)) set.seed(isim)

    ## Generate data + path
    maxsteps = 5
    y0 = fun(lev,n)
    D = genlassoinf::makeDmat(n,type='tf',ord=0)

    ## Obtain stoptime
    f0    = dualpathSvd2(y0, D=D, maxsteps, approx=T)
    mm = get.modelinfo(obj=f0, consec=2, sigma=sigma)
    bic = mm$ic
    stopped = mm$stopped
    while(!stopped){
        f0    = dualpathSvd2(y0, D=D, maxsteps, approx=T)
        mm = get.modelinfo(obj=f0, consec=2, sigma=sigma)
        bic = mm$ic
        stopped = mm$stopped
        maxsteps = maxsteps + 5
    }

    ## 2-rise rule
    f1 = stop_path(f0, sigma=sigma, stoprule="bic")
    stoppedmodel = f1$stoppedmodel
    stoppedmodel = stoppedmodel[!(is.na(stoppedmodel))]
    numcp.tworise = length(stoppedmodel)

    ## Global minimization rule.
    numcp.global = which.min(bic)-1
    numcp.global = pmin(numcp.global, n-consec-1)

    return(list(numcp.tworise=numcp.tworise,
                numcp.global= numcp.global))
}


## Run the onejump example
nlev = 6
levs = seq(from=0.01,to=3,length=nlev)
nsim = 2000
stoptimelist = mclapply(1:length(levs), function(ilev){
    print(ilev)
    lev = levs[ilev]
    results <- replicate(nsim, {onesim(lev,n,onejump)})
    return(results)
},mc.cores=3)

## Save this data
settings = list(nsim=nsim,levs=levs,onesim=onesim,levs=levs)
save(list=c("settings", "stoptimelist"), file="../output/ic-performance.Rdata")

## Percent of times the right thing was captured
getstopprop = function(stoptimemat, njumps){
    prop1 = unlist(stoptimemat[1,])
    prop2 = unlist(stoptimemat[2,])
    prop1 = sum(unlist(stoptimemat[1,]) %in% njumps)/nsim
    prop2 = sum(unlist(stoptimemat[2,]) %in% njumps)/nsim
    return(c(global = prop1, tworise = prop2))
}
getstopprop12 = function(stoptimemat){getstopprop(stoptimemat, c(1,2))}
stopprop12 = t(sapply(stoptimelist, getstopprop12))

## Make plot
w=5
h=5
pdf(file="../output/ic-performance.pdf", width=w, height=h)
lty=c(1,2)
legend = c("Two-rise", "Global")
col=c(1,1)
plot(NA, ylim=c(0,1), xlim=c(0,max(levs)),
     axes= FALSE, ylab="lev", xlab= "Proportion")
axis(1)
axis(2)
matlines(y=(stopprop12), x=levs,  type='l', lty=lty, col=col)
legend("bottomright", lty=lty, legend = legend)
graphics.off()

## Bands; no need to plot, since they are so small.
upmat = downmat = matrix(NA,ncol=ncol(stopprop12),nrow = nrow(stopprop12))
for(ii in 1:nrow(stopprop12)){
    for(jj in 1:ncol(stopprop12)){
        prop = stopprop12[ii,jj]
        band = 1/nsim*prop*(1-prop)
        upmat[ii,jj] = stopprop12[ii,jj] + band
        downmat[ii,jj] = stopprop12[ii,jj] - band
    }
}
