## Synopsis: Calculate conditional power after one or two steps, in an
## alternating jump example, for the bottom row of figure 8. Output directly
## used by ./1d-alternjump.R
outputdir = "../output"

## Simulation settings
lev1= 0
lev2 = 2
sigma = 1
maxcount = 1000
alpha = 0.05
numsteps = 2
isim = 0

## icount and jcount are for counting the number of first&second step
## conditional cases encountered
icount = jcount = 0
p1spike = p1segment = p21spike = p21segment = rep(NA, maxcount)
ds1segment = ds1spike = ds21segment = ds21spike = list()
cis1.segment = cis21.segment = cis1.spike = cis21.spike = matrix(NA,nrow=maxcount,ncol=2)
cat("Simulations running with lev2=", lev2, fill=TRUE)
while(icount < maxcount & jcount < maxcount){

    printprogress(icount, maxcount)
    isim = isim + 1

    ## Generate data and obtain path + polyhedron
    y0 <- alternjump.y(lev1=lev1, lev2=lev2, sigma=sigma)
    beta0 <- alternjump.y(returnbeta=T, lev1=lev1, lev2=lev2, sigma=sigma)
    f1 <- dualpathSvd2(y0, dual1d_Dmat(length(y0)), maxsteps=1,verbose=FALSE,approx=TRUE)
    f2 <- dualpathSvd2(y0, dual1d_Dmat(length(y0)), maxsteps=2,verbose=FALSE,approx=TRUE)


    ## Form conditional p-values for first step, from the first step
    if(f1$pathobj$B[1] %in% c(20,40)){

     icount=icount+1

        ## Form contrasts
        G = f1$Gobj.naive$G
        u = f1$Gobj.naive$u
        d1spike <- getdvec(obj=f1,y=y0,k=1,type="spike",scale="segmentmean")
        d1segment <- getdvec(obj=f1,y=y0,k=1,type="segment",scale="segmentmean")

        ## Form p-values & one-sided confidence intervals
        p1spike[icount] <- poly.pval(y=y0, G=G, v=d1spike, u=u, sigma=sigma)$pv
        p1segment[icount] <- poly.pval(y=y0, G=G, v=d1segment, u=u, sigma=sigma)$pv
        cis1.segment[icount,] <- confidence_interval(y0,list(gamma=G,u=u),
                                                  contrast=d1segment, sigma=sigma,
                                                  alpha=alpha,
                                                  alternative="one.sided",
                                                  fac=10)
        cis1.spike[icount,] <- confidence_interval(y0,list(gamma=G,u=u),
                                                  contrast=d1spike, sigma=sigma,
                                                  alpha=alpha,
                                                  alternative="one.sided",
                                                  fac=10)
        ds1segment[[icount]] = d1segment
        ds1spike[[icount]] = d1spike
    }

    ## Form conditional p-values for first jump, from the second step
    if( all(f2$pathobj$B[1:2] %in%  c(20,40))){
        jcount=jcount+1

        ## Form contrasts
        G = f2$Gobj.naive$G
        u = f2$Gobj.naive$u
        d21spike <- getdvec(obj=f2,y=y0,k=1,type="spike")
        d21segment <- getdvec(obj=f2,y=y0,k=1,klater=2,type="segment")

        ## Form p-values & one-sided confidence intervals
        p21spike[icount] <- poly.pval(y=y0, G=G, v=d21spike, u=u, sigma=sigma)$pv
        p21segment[icount] <- poly.pval(y=y0, G=G, v=d21segment, u=u, sigma=sigma)$pv
        cis21.segment[icount,] = confidence_interval(y0,list(gamma=G,u=u),
                                                  contrast=d21segment,
                                                  sigma=sigma, alpha=alpha,
                                                  alternative="one.sided",
                                                  fac=10)
        cis21.spike[icount,] = confidence_interval(y0,list(gamma=G,u=u),
                                                  contrast=d21spike,
                                                  sigma=sigma, alpha=alpha,
                                                  alternative="one.sided",
                                                  fac=10)
        ds21segment[[icount]] = d21segment
        ds21spike[[icount]] = d21spike
    }
}

print(icount)
## Trim things
results = list(cis1.segment = cis1.segment[!is.na(p1segment),],
               cis21.segment = cis21.segment[!is.na(p21segment),],
               cis1.spike = cis1.spike[!is.na(p1segment),],
               cis21.spike = cis21.spike[!is.na(p21segment),],
                  p1spike = p1spike[!is.na(p1spike)],
                  p21spike = p21spike[!is.na(p21spike)],
                  p1segment = p1segment[!is.na(p1segment)],
                  p21segment = p1segment[!is.na(p21segment)], icount = icount-1,
                  jcount = jcount-1,
                  propcorrect.step1 = icount/isim,
                  propcorrect.step2 = jcount/isim,
                  ds1spike=ds1spike[!is.na(p1spike)],
                  ds1segment=ds1segment[!is.na(p1spike)],
                  ds21spike=ds21spike[!is.na(p21spike)],
                  ds21segment=ds21segment[!is.na(p21segment)])
print(results$icount)

filename = paste0("updown-example-lev",lev2,".Rdata")
print(filename)
save(results, file=file.path(outputdir,filename))
