## Make sure you're working from [dropboxfolder]/code
source('selectinf/selectiveInference/R/funs.inf.R')
source('funs.R')
source('testfuns.R')
source('dualPathSvd2.R')
library(genlasso)
library(igraph)
verbose = F
outputdir = "output"
codedir = "."

## Generate simulation results
nsim = 100000
ngrain=5+1
lev1 = 0
lev2 = 0
sigma=15/15
nn = 10
Dmat = graph2D_Dmat(nn^2)
maxsteps = 100
consec = 1
mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat), mode="undirected")
igrain=6 # the sixth grain is the null
ngrain=6
icount = 0
nhit.enough = 3000
pvals = matrix(NA,nrow=nsim,ncol=ngrain)
cat("Running simulation with nn x nn image, nn= ", nn, "for lev2=",lev2, fill=T)
start.time = Sys.time()
for(isim in 1:nsim){
    if(icount >nhit.enough) break
    curr.time = Sys.time()
            cat("\r", icount, "out of", nhit.enough, "and", round(as.numeric(curr.time-start.time))/60, "minutes have passed.")

        # Make a noisy nn x nn image
            A = matrix(lev1,ncol=nn, nrow = nn)
            A[1:(nn/2),1:(nn/2)] = lev2
            A.noisy = A + rnorm(nn^2,mean=0,sd=sigma)

        # generate data
            beta0 = do.call(c,lapply(1:nrow(A),function(irow)A[irow,]))
            y0 = do.call(c,lapply(1:nrow(A.noisy),function(irow)A.noisy[irow,]))
            # optionally, plot!
            image(A.noisy)

        # fit graph fused lasso using the matrix given by D only (initially go a generous number of steps)
            f0 = dualpathSvd2(y0,Dmat,verbose=FALSE, maxsteps = maxsteps)
            maxsteps = length(f0$action[f0$action!=0])

        # Get BIC scores
            mm = get.modelinfo(f0,y0,sigma,maxsteps,Dmat, stoprule = 'bic')
            bic = mm$ic
            stop.time = which.rise(bic,consec)-1

        # Get polyhedron (with BIC)
            if(stop.time+consec < maxsteps){
            Gobj.new.with.stoptime = getGammat.with.stoprule(obj=f0,y=y0,
                                            condition.step = stop.time+consec,
                                            stoprule = "bic", sigma=sigma,
                                            consec=consec, maxsteps=maxsteps, type = "graph", D = Dmat)
            } else {
            print('bic rule hasnt stopped!')
            next
            }
            G = Gobj.new.with.stoptime$Gammat
            u = Gobj.new.with.stoptime$u

            if(!polyhedron.checks.out(y0,G,u)){
                print("polyhedron is problematic")
                next
            }

        # Calculate p-value and store it
        steps.to.test = mm$knots.primal[mm$knots.primal<=stop.time]
        steps.to.test = steps.to.test[steps.to.test!=1]
        states = get.states(f0$action)
        final.model = states[[stop.time+1]]
        final.model.signs = f0$ss[[stop.time+1]]

        for( mystep  in steps.to.test ){

            # Get the residual subspace basis
            v = mm$resids[,mystep]

            # Adjust it by the sign
            test.knot = final.model[mystep-1]
            slopedir = sign(Dmat[test.knot,] %*% v)
            test.knot.sign = final.model.signs[mystep-1]
            v = v * test.knot.sign * sign(slopedir-.5) 

            ## if(any(is.na(v))) print(mm)
            icount = icount + 1
            pvals[icount,igrain] = pval.fl1d(y0, G, v, sigma, u = u)
        }
    save(ngrain, nsim, Dmat, maxsteps, mygraph, pvals, nhit.enough,icount, igrain,
     file = file.path(outputdir, paste0("paper-2dgraph-null.Rdata")))
}
cat(fill=T)
