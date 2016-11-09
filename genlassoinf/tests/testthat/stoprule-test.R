context("Test polyhdron forming")

sigma = .1
n = 100
maxsteps = 5
lev1 = 0
lev2 = 3
consec=2
D = dual1d_Dmat(length(y0))

# Generate data + path
beta0 = rep(c(lev1,lev2),each=n/2)
y0    = beta0 + rnorm(n, 0,sigma)
f0    = dualpathSvd2(y0, D=D, maxsteps, approx=T)
mm = get.modelinfo(f0, y0, sigma, maxsteps, D=D)
bic = mm$ic  
stoptime.bic = which.rise(bic,consec) - 1
stoptime.bic = pmin(stoptime.bic, n-consec-1)
locs.bic = f0$pathobj$B[1:stoptime.bic]
Gobj    = getGammat.with.stoprule(obj=f0,y=y0,
                                  condition.step=stoptime.bic+consec, type='tf',
                                  stoprule='bic',sigma=sigma,consec=consec,
                                  maxsteps=maxsteps, D=D)

test_that("After appending stopping rule, Gamma matrix and u vector have appropriate dimensions", {
    expect_equal(nrow(Gobj$Gammat) , length(Gobj$u))
}

