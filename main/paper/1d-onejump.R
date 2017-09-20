## Synopsis: For one-jump example in figure 8, generate data and plots and
## confidence intervals

outputdir = "../output"

## Generate p-values and simulation quantities
nsim = 1000
n = 60
sigma = 1
lev1= 0
ngrain = 3
lev2list = seq(from=0,to=2, length(ngrain))#c(0,1,2)
numsteps=1
spike = onejump.naive.sim(testtype = "spike", nsim=nsim,sigma=sigma,lev1=lev1,
                          lev2list=lev2list,numsteps=numsteps,verbose=T, loctype="exact")
segment = onejump.naive.sim(testtype = "segment", nsim=nsim,sigma=sigma,lev1=lev1,
                            lev2list=lev2list,numsteps=numsteps,verbose=T, loctype="exact")
save(file.path(outputdir,"onejump-example.Rdata"))


## Plot settings
load(file.path(outputdir,"onejump-example.Rdata"))
set.seed(0)
sigma = 1
y0    = onejump.y(returnbeta=F,lev1=0,lev2=2,sigma=sigma,n=60)
beta0 = lapply(c(0:2), function(lev2) onejump.y(returnbeta=T,lev1=0,lev2=lev2,sigma=sigma,n=60))
beta0.middle = onejump.y(returnbeta=T,lev1=0,lev2=1,sigma=sigma,n=60)
beta0.top = onejump.y(returnbeta=T,lev1=0,lev2=2,sigma=sigma,n=60)
x.contrasts = c(1:60)
v.spike = c(rep(NA,29),.5*c(-1,+1)-2, rep(NA,29))
v.segment = c(.3*rep(-0.7,30)-2 , .3*rep(+0.7,30)-2 )

xlab = "Location"
w = 5; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-3,5)
mar = c(4.5,4.5,0.5,0.5)
xlim = c(0,70)

xticks = c(0,2,4,6)*10
let = c("A","B")
ltys.sig = c(2,2,1)
lwd.sig = 2
pch.dat = 16
pcol.dat = "grey50"
pch.contrast = 17
lty.contrast = 2
lcol.sig = 'red'
pcol.spike=3
pcol.segment=4
pcols.delta =   pcols.oneoff = RColorBrewer::brewer.pal(n=3,name="Set2")
pch.spike = 15
pch.segment = 17
cex.contrast = 1.2

## Plot single example of data and contrast (top leftmost figure)
pdf(file.path(outputdir,"onejump-example-data-and-contrast.pdf"), width=5,height=5)

par(mar=c(4.1,3.1,3.1,1.1))
 plot(y0, ylim = ylim,axes=F, xlim=xlim, xlab = xlab, ylab = "", pch=pch, col=pcol);
 axis(1, at = xticks, labels = xticks); axis(2)
 for(ii in 1:3) lines(beta0[[ii]],col="red",lty=ltys.sig[ii], lwd=lwd.sig)
 for(ii in 0:2) text(x=65,y=ii, label = bquote(delta==.(ii)))
 points(v.spike~x.contrasts, pch = pch.spike, col = 3)
 points(v.segment~x.contrasts, pch = pch.segment, col = 4)
 abline(h = mean(v.segment,na.rm=T), col = 'lightgrey')
 legend("topleft", pch=c(pch.dat,NA,pch.spike,pch.segment),
        lty=c(NA,1,NA,NA), lwd=c(NA,2,NA,NA),
        col = c(pcol.dat, lcol.sig, pcol.spike,pcol.segment),
        pt.cex = c(cex.contrast, NA, cex.contrast, cex.contrast),
        legend=c("Data", "Mean","Spike contrast", "Segment contrast"))
title(main=expression("Data example"))
graphics.off()

## Create QQ plot of correct location p-values (top middle and right figure)
contrast.type = c("spike", "segment")
dat = list(spike=spike, segment=segment)
for(jj in 1:2){
    pdf(file.path(outputdir,paste0("onejump-example-qqplot-",tolower(contrast.type[jj]),".pdf")),
        width=5,height=5)
    mydat = dat[[contrast.type[jj]]]
    unif.p = seq(from=0,to=1,length=nsim)
    for(ii in 1:3){
        if(ii!=1) par(new=T)
        a = qqplot(x=unif.p, y=mydat$pvals.correctlist[[ii]], plot.it=FALSE)
        if(ii==1){
            plot(x=a$y, y=a$x, axes=FALSE,xlab="",ylab="", col = pcols.delta[ii], pch=16)
        } else {
            points(x=a$y, y=a$x, col = pcols.delta[ii], pch=16)
        }
    }
    axis(2);axis(1)
    abline(0,1)
    mtext("Expected",1,padj=4)
    mtext("Observed",2,padj=-4)
    title(main = bquote(.(contrast.type[jj])~test~p-values))
    if(jj==1){
      legend("bottomright", col = pcols.delta,
             lty = 1, lwd = 5,
             legend = sapply(c(bquote(delta==0),
                               bquote(delta==1),
                               bquote(delta==2)), as.expression) )
    }
    graphics.off()
}


## Create Confidence interval coverage tables
contrast.type = c("spike", "segment")
dat = list(spike=spike, segment=segment)
twocoverages = lapply(1:2, function(jj){
    mydat = dat[[contrast.type[jj]]]
    coverages = sapply(1:ngrain, function(igrain){
        ci.list = (dat[[jj]])$cis.correctlist[[igrain]]
        ds = dat[[jj]]$dlist[[igrain]]
        lev2 = lev2list[[igrain]]
        ## truejumpsize = sum(d * onejump(lev=lev2list[[igrain]],n=n))
        ## coverage = sum(sapply(ci.list, function(myci) return(myci[1] < truejumpsize)))/nsim
        covered = Map(function(myci, myd){
            truejumpsize = sum(myd * onejump(lev=lev2,n=n))
            return(myci[1] < truejumpsize)
        }, ci.list,ds)
        coverage = sum(unlist(covered))/nsim
        return(coverage)
    })
    coverages = rbind(lev2list,coverages)
    rownames(coverages) = c("delta", "coverages")
    return(coverages)
})
names(twocoverages) = contrast.type
xtable::xtable(twocoverages[["spike"]])
xtable::xtable(twocoverages[["segment"]])




## Get average (actually, median) lower bounds of one-sided confidence intervals
twolowerbounds = lapply(1:2, function(jj){
    mydat = dat[[contrast.type[jj]]]
    av.lowers = sapply(1:ngrain, function(igrain){
        ci.list = (dat[[jj]])$cis.correctlist[[igrain]]
        all.lowers = sapply(ci.list, function(myci){
            myci[1]
        })
        mdn = median(all.lowers)
        ## avg = median(all.lowers)
        return(mdn)
    })
    mdn.lowers = rbind(lev2list,round(av.lowers,3))
    rownames(mdn.lowers) = c("delta", "mdn.lowers")
    return(mdn.lowers)
})


