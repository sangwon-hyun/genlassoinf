context("Testing whether path class works well, on linear trend filtering example.")

## Simulation settings
beta1 = rep(0,20)
fac = .1
beta2.low = .5*fac*seq(from=21,to=40) - 10*fac
beta0 = c(beta1,beta2.low)
n = length(beta0)

## Run simulation
consec = 2
sigma = 1
tf.order = 1
D = makeDmat(n,order=tf.order) 
maxsteps = 20

one_check = function(){

    ## Generate Data and path
    y0 = beta0 + rnorm(length(beta0),0,sigma)
    f0 = dualpathSvd2(y0, D, maxsteps = maxsteps,verbose=F)
    
                                        # Collect Gammat at stop time
    bic   = get.modelinfo(f0,y0,sigma,maxsteps,D=D, stoprule = 'bic')$ic
    stop.time = which.rise(bic,consec=consec) - 1
    stop.time = pmin(stop.time,n-consec-1)
    
    if(!(stop.time+consec < maxsteps)){
        print('bic rule hasnt stopped!')
        next
    }
    
    Gobj.new.with.stoptime = getGammat.with.stoprule(obj=f0,y=y0, condition.step
                                                     = stop.time+consec, stoprule = "bic",
                                                     sigma=sigma, type='tf', consec=consec,
                                                     maxsteps=maxsteps,D=D)
    G = Gobj.new.with.stoptime$Gammat
    u = Gobj.new.with.stoptime$u

    return(polyhedron.checks.out(y0,G,u))
}

check.results = replicate(10, one_check())

expect_true(all(check.results))
