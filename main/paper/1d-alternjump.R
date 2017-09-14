## Synopsis: For up-then-down-jump example in figure 8, generate data and plots
## and /conditional/ confidence intervals, given that the right location was
## selected.

outputdir = "../output"

### Up-then-down signal Example ######

## Generate p-values and simulation quantities
#  # Generate stuff
#  lev1  = 0
#  sigma = 1
#  nsim = 10000
#

## Plot settings
load(file=file.path(outputdir,"updown-example.Rdata"))
## oldoutputdir = "/home/justin/Dropbox/research/genlassoinf/code/output"
## oldoutputdir = "/media/shyun/Bridge/Dropbox/research/genlassoinf/code/output"
## load(file=file.path(oldoutputdir,"updown-example.Rdata"))


w=5
h=5
xlab = "Location"
w = 5; h = 5
pch = 16; lwd = 2
pcol = "gray50"
ylim = c(-6,5)
mar = c(4.5,4.5,0.5,0.5)
xlim = c(0,70)
xticks = c(0,2,4,6)*10
let = c("A","B")
ltys.sig = 1
lwd.sig = 2
pch.dat = 16
pcol.dat = "grey50"
pch.contrast = 17
lty.contrast = 2
lcol.sig = 'red'
pcol.spike=3
pch.qq=16
pcol.segment=4
pcols.delta = RColorBrewer::brewer.pal(n=3,name="Set2")
pch.segment = 17
pch.spike = 15
cex.contrast = 1.2

## Plot single example of data and contrast (bottom leftmost figure)
pdf("output/updown-example-data-and-contrast.pdf", width=w,height=h)
par(mar=c(4.1,3.1,3.1,1.1))
set.seed(0)
sigma = 1
y0    = alternjump.y(returnbeta=F,lev1=0,lev2=2,sigma=sigma,n=60)
beta0 = alternjump.y(returnbeta=T,lev1=0,lev2=2,sigma=sigma,n=60)

x.contrasts = c(1:60)
v.spike = c(rep(NA,29),.5*c(-1,+1)-2, rep(NA,29))
v.segment = c(.3*rep(-0.7,30)-2 , .3*rep(+0.7,30)-2 )

v1.spike = c(rep(NA,19),.6*c(-1,+1)-2, rep(NA,39))-1
v1.segment = c(.3*rep(-0.7,20)-2 , .3*rep(+0.7,20)-2 , rep(NA,20))-1

v2.spike = c(rep(NA,19),.6*c(-1,+1)-2, rep(NA,39))-3
v2.segment = c(.3*rep(-0.7,20)-2 , .3*rep(+0.35,40)-2)-3


plot(y0, ylim = ylim,axes=F, xlim=xlim, xlab = xlab, ylab = "", pch=pch, col=pcol);
axis(1, at = xticks, labels = xticks); axis(2)
lines(beta0,col="red",lty=ltys.sig, lwd=lwd.sig)
ii=2; text(x=45,y=ii, label = bquote(delta==.(ii)))
points(v1.segment~x.contrasts, pch = pch.segment, col = 4)
points(v1.spike~x.contrasts, pch = pch.spike, col = 3)
points(v2.segment~x.contrasts, pch = pch.segment, col = 4)
points(v2.spike~x.contrasts, pch = pch.spike, col = 3)

lines(x = c(0,61), y = rep(mean(v1.segment,na.rm=T),2), col = 'lightgrey' )
lines(x = c(0,61), y = rep(mean(v2.segment,na.rm=T),2), col = 'lightgrey' )
ii=1; text(x=67,y=mean(v1.segment,na.rm=T), label = bquote(Step~.(ii)))
ii=2; text(x=67,y=mean(v2.segment,na.rm=T), label = bquote(Step~.(ii)))

title(main=expression("Data example"))
graphics.off()

## QQ plot of correct location p-values (bottom, middle and right plot)
pspikes = list(p1spike, p21spike)
psegments = list(p1segment, p21segment)
pvals.list = list(pspikes, psegments)

w=5
h=5
contrast.type = c("Spike", "Segment")
for(jj in 1:2){
    ## pdf(paste0("output/updown-example-qqplot-",tolower(contrast.type[jj]),".pdf"), width=5,height=5)
    pdf(file.path(outputdir, paste0("updown-example-qqplot-",tolower(contrast.type[jj]),".pdf")), width=w,height=h)
    unif.p = seq(from=0, to=1, length=nsim)
    for(ii in 1:2){
        pvals = pvals.list[[jj]]
        if(ii!=1) par(new=T)
        a = qqplot(y=pvals[[ii]], x=unif.p, plot.it=FALSE)
        myfun = (if(ii==1) plot else points)
        suppressWarnings(
            myfun(x=a$y, y=a$x, axes=F, xlab="", ylab="", col = pcols.delta[ii], pch=pch.qq)
        )
        abline(0,1,col='lightgrey')
    }
    axis(2);axis(1)
    mtext("Observed",2,padj=-4)
    mtext("Expected",1,padj=4)
    title(main = bquote(.(contrast.type[jj])~test~p-values))
    if(jj==1){
        legend("bottomright", col = pcols.delta,
               lty = 1, lwd = 5,
               legend = sapply(c(bquote(Step~1),
                                 bquote(Step~2)), as.expression) )
    }
    graphics.off()
}



## Create Confidence interval coverage tables (change appropriately)
contrast.type = c("spike", "segment")
dat = list(spike=spike, segment=segment)

mn = alternjump.y(lev1=lev1, lev2=lev2, sigma=sigma,returnbeta=TRUE)
coverage1 = sum(sapply(1:icount, function(ii){
    myv = ds1[[ii]]
    myci = cis1.segment[ii,]
    truejumpsize = sum(myv*mn)
    return(myci[1]<truejumpsize)
}))/icount

coverage21 = sum(sapply(1:icount, function(ii){
    myv = ds1[[ii]]
    myci = cis1.segment[ii,]
    truejumpsize = sum(myv*mn)
    return(myci[1]<truejumpsize)
}))/icount


sum(apply(cis21.segment, 1, function(myci) return(myci[1] < truejumpsize)))/nrow(cis21.segment)

sum(apply(cis1.segment, 1, function(myci) return(myci[1] < truejumpsize)))/nrow(cis1.segment)


twocoverages = lapply(1:2, function(jj){
    mydat = dat[[contrast.type[jj]]]
    coverages = sapply(1:ngrain, function(igrain){
        ci.list = (dat[[jj]])$cis.correctlist[[igrain]]
        d = dat[[jj]]$d[[igrain]]
        truejumpsize = 2
        coverage = sum(sapply(ci.list, function(myci) return(myci[1] < truejumpsize)))/nsim
        return(coverage)
    })
    coverages = rbind(lev2list,coverages)
    rownames(coverages) = c("delta", "coverages")
    return(coverages)
})
names(twocoverages) = contrast.type
xtable::xtable(twocoverages[["spike"]])
xtable::xtable(twocoverages[["segment"]])
