## Generate data
sigma = .1
n = 100
maxsteps = 5
lev1 = 0
lev2 = 3
consec=2
D = dual1d_Dmat(length(y0))

## Generate data + path
beta0 = rep(c(lev1,lev2),each=n/2)
y0    = beta0 + rnorm(n, 0,sigma)
f0    = dualpathSvd2(y0, D=D, maxsteps, approx=T)

## Get naive poyhedron (fixed stop time of 1)
Gobj    = getGammat.naive(obj = f0, y = y0, condition.step = stoptime+consec)
G = Gobj$Gammat
u = Gobj$u

## Form contrast and segment test p-value
states = get.states(f0$action)
stoptime = 1
final.model = states[[stoptime+1]]
for(ii in 1:length(final.model)){
    this.sign = f0$pathobj$s[which(f0$pathobj$B == final.model[ii])]
    my.v.lrt = make.v.tf.fp(test.knot = final.model[ii],
                            adj.knot  = final.model,
                            test.knot.sign = this.sign,
                            D=D)
    pval = poly.pval(y=y0, G=G, u=u, v=v, sigma=sigma)$pv
    cat("After fixed size model was selected, the segment test pvalue at location",
        final.model[ii], "is", pval,fill=TRUE)
}

## Get stop-time-incorporated polyhedron.
mm = get.modelinfo(f0, y0, sigma, maxsteps, D=D)
bic = mm$ic  
stoptime = which.rise(bic,consec) - 1
stoptime = pmin(stoptime, n-consec-1)
Gobj    = getGammat.with.stoprule(obj = f0, y = y0,
                                  condition.step = stoptime+consec, type ='tf',
                                  stoprule = "bic", sigma = sigma,
                                  consec = consec, maxsteps = maxsteps, D = D)
G = Gobj$Gammat
u = Gobj$u

## Form contrast and segment test p-values
final.model = states[[stoptime+1]]
for(ii in 1:length(final.model)){
    this.sign = f0$pathobj$s[which(f0$pathobj$B == final.model[ii])]
    my.v.lrt = make.v.tf.fp(test.knot = final.model[ii],
                            adj.knot  = final.model,
                            test.knot.sign = this.sign,
                            D=D)
    pval = poly.pval(y=y0, G=G, u=u, v=v, sigma=sigma)$pv
    cat("After model selected by BIC-rise stop rule, the segment test pvalue at location",
        final.model[ii], "is", pval,fill=TRUE)
}
