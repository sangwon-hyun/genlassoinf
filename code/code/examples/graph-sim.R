## Make sure you're working from [dropboxfolder]/code
#  workingdir = '~/../../media/shyun/Bridge/Dropbox/CMU/courses(CURRENT)/genlassoinf/code'
#  setwd(workingdir)
  source("settings.R")
  source('selectinf/selectiveInference/R/funs.inf.R')
  source('funs.R')
  source('examples/testfuns.R')
  source('dualPathSvd2.R')
  library(genlasso)
  library(igraph)
  verbose = F

# simulation settings
#  nsim = 500
#  ngrain = 5
#  sigmamin = 2
#  sigmamax = 4
#  sigmalist = seq(from = sigmamin, to = sigmamax, length = ngrain)
  nsim = 1
  #sigmalist = c(0.1, 0.5, 1)#seq(from = 0.1, to = 2, length = ngrain)
  sigma=1
  ngrain = length(sigmalist)
  maxsteps = 30
  
  chunks = list(1,2,3,4)#list(lev2list[1:3],lev2list[4:5])
  nchunk = 4
  
 
# generate some basic things
  # make adjacency matrix
  set.seed(100)
  clustersize = 20
  nclust = 3
  corr = 0.01
  block.cor.mat = matrix(c(1,corr,corr,
                           corr,1,corr,
                           corr,corr,1), nrow=3)
  mygraph = sample_sbm(clustersize*3, pref.matrix = block.cor.mat, block.sizes = rep(clustersize,3))
  
# make D matrix from the adjacency matrix
  adjmat = get.adjacency(mygraph)
  D = getDmat.from.adjmat(adjmat)
  
# get ready to store things
  pmat = array(dim = c(nsim,ngrain,6))
  stoptimemat = array(dim = c(nsim, ngrain,6))
  pairindlist = list(1:2,2:3,c(3,1),c(1,23), c(2,13), c(3,12))
  pairindlist = lapply(pairindlist, function(mypair){if(mypair[2]>10) return(c(mypair[1], paste(floor(mypair[2]/10), " and ", mypair[2]%%10  ))   ) else return(mypair) })
  consec=1  
  
  

# main loop
for(igrain in 1:ngrain){

  cat(igrain, "out of", ngrain, fill=T)
  sigma = sigmalist[igrain]
  maxsteps = 300#min(100 + igrain*5, 120)

  for(isim in 1:nsim){

    cat("\r", isim, "out of", nsim)
 # Make a noisy nn x nn image

      # generate data
        beta0 = c(rep(0,clustersize), rep(1,clustersize), rep(3,clustersize))
        y0 = beta0 + rnorm(clustersize*3,0,sigma)

      # fit graph fused lasso using the matrix given by D only (initially go a generous number of steps)
        f0 = dualpathSvd2(y0,D,verbose=T, maxsteps = maxsteps)
        
      # Get BIC scores
        mm = get.modelinfo(f0,y0,sigma,maxsteps,D, stoprule = 'bic')
        bic = mm$ic
        stop.time = which.rise(bic,consec)-1
        print(stop.time)
        if(stop.time >= maxsteps){
          print(paste("didn't stop by step", maxsteps));
          next
        }

      # Get polyhedron (with BIC)
        if(stop.time+consec < maxsteps){
        Gobj.new.with.stoptime = getGammat.with.stoprule(obj=f0,y=y0,
                                           condition.step = stop.time+consec,
                                           stoprule = "bic", sigma=sigma,
                                           consec=consec, maxsteps=maxsteps, type = "graph", D = D)
        } else {
          print('bic rule hasnt stopped!')
          next
        }
        G = Gobj.new.with.stoptime$Gammat
        u = Gobj.new.with.stoptime$u
        
        #G%*%y0 - u
        
        
#        ## Temporary check if BIC rule is correct
#        for(iisim in 1:10){
#        print(iisim)
#          iii = which(u!=0)
#          y1 = y0 + rnorm(n,0,0.000001)
#          f1 = dualpathSvd2(y1,D,verbose=FALSE, maxsteps = maxsteps)
#          mm1 = get.modelinfo(f1,y1, sigma, maxsteps, D)
#          if(which.rise(mm1$ic,consec) - 1 == 19){
#            print(all((G[iii,]%*%y1 - u[iii])>0))
#          }
#        }

        # Check correctness of polyhedron:
        polyhedron.checks.out = function(y0,G,u){
          all.nonneg = all((G%*%y0 >= u))
          if(all.nonneg){ 
            return(all.nonneg) 
          } else {
            #print(tail(sort(G%*%y0 - u, decreasing=T)))
#            print("End of G * y0 - u")
#            print(tail(G%*%y0 - u))
#            stop("failed Gy >= u test")
            print("failed Gy >= u test")
            return(all.nonneg)
          }
        }        
        #stopifnot(polyhedron.checks.out(y0,G.bic,u.bic))
        if(!polyhedron.checks.out(y0,G,u)){
            print("polyhedron is problematic")
            next
        }

        

     # Calculate p-value and store it
       steps.to.test = mm$knots.primal[mm$knots.primal<=stop.time & mm$actiondirs.primal ==1]
       steps.to.test = steps.to.test[steps.to.test!=1]
       steps.to.test = steps.to.test[!is.na(steps.to.test)]
       for( mystep  in steps.to.test ){
         
         v = mm$resids[,mystep]
         if(any(is.na(v))) print(mm)
        # if( v %*% y0 < 0) v = -v

           # Helper function to see if the segment contrast is correct, in the graph case
           contrast.is.right=function(v){
             factors = factor(beta0)
             unique.factors = levels(factor(beta0))
             ind.1 = which(factors == unique.factors[1])
             ind.2 = which(factors == unique.factors[2])
             ind.3 = which(factors == unique.factors[3])
             v = unname(v)
             v = round(v,5)
             
             g1.g2 = (all.equal(v[ind.1],rep(v[ind.1[1]], length(ind.1)))==T) & (all.equal(v[ind.2],rep(v[ind.2[1]], length(ind.2)))==T) & all(v[ind.1] != v[ind.2]) & (all.equal(v[!(1:length(v) %in% c(ind.1,ind.2))], rep(0,length(v)-length(ind.1)-length(ind.2)))==T)
             g2.g3 = (all.equal(v[ind.2],rep(v[ind.2[1]], length(ind.2)))==T) & (all.equal(v[ind.3],rep(v[ind.3[1]], length(ind.3)))==T) & all(v[ind.2] != v[ind.3]) & (all.equal(v[!(1:length(v) %in% c(ind.2,ind.3))], rep(0,length(v)-length(ind.2)-length(ind.3)))==T)
             g3.g1 = (all.equal(v[ind.3],rep(v[ind.3[1]], length(ind.3)))==T) & (all.equal(v[ind.1],rep(v[ind.1[1]], length(ind.1)))==T) & all(v[ind.1] != v[ind.3]) & (all.equal(v[!(1:length(v) %in% c(ind.1,ind.3))], rep(0,length(v)-length(ind.1)-length(ind.3)))==T)

             g1.g23 = (all.equal(v[ind.1],rep(v[ind.1[1]], length(ind.1)))==T) & (all.equal(v[c(ind.2,ind.3)],rep(v[ind.2[1]], length(c(ind.2,ind.3))))==T)
             g2.g13 = (all.equal(v[ind.2],rep(v[ind.2[1]], length(ind.2)))==T) & (all.equal(v[c(ind.1,ind.3)],rep(v[ind.1[1]], length(c(ind.1,ind.3))))==T) 
             g3.g12 = (all.equal(v[ind.3],rep(v[ind.3[1]], length(ind.3)))==T) & (all.equal(v[c(ind.1,ind.2)],rep(v[ind.1[1]], length(c(ind.1,ind.2))))==T)
                        
             g.all = c(g1.g2,g2.g3,g3.g1,g1.g23,g2.g13,g3.g12)           
             return(g.all)
           }

           # Wrapper for contrast.is.right()
           which.comparison = function(v){which(contrast.is.right(v))}

           # Conduct comparison of contrast is right
           if(any(contrast.is.right(v))){
             pmat[isim,igrain, which.comparison(v)] = pval.fl1d(y0, G, v, sigma, u = u)
             stoptimemat[isim,igrain,which.comparison(v)] = stop.time
             if(is.nan(pmat[isim,igrain,which.comparison(v)])) stop("nan pvalue occurred")
         }
      }
   }
   cat(fill=T)
   save(ngrain,nsim,sigmalist,clustersize,nclust,adjmat,
        D,maxsteps,mygraph,corr,block.cor.mat,pmat,
        stoptimemat, pairindlist, igrain, 
        file = file.path(outputdir, paste0("sbm-20each-chunk-",ichunk,".Rdata")))
 }



 
  ## See what is happening with the breaks
#  pdf("~/Desktop/see-sbm-breaks.pdf",width=20,height=20)
#  par(mfrow=c(5,5))
#    for(step in 1:min(stop.time,25)){
#      print(step)
#      Dmat.curr = Dmat[-getB.from.actions(f0$action[1:step]),]
#      mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat.curr), mode="undirected")
#      plot(mygraph)
#      title(main=step)
#      }
#  dev.off()


       
### Also do some basic checks: is the polyhedron without error?
#for(jj in 1:5){
#  sigma=seq(from=0.1,to=10,length=5)[jj]
#  beta0 = c(rep(0,clustersize), rep(1,clustersize), rep(3,clustersize))
#  y0 = beta0 + rnorm(clustersize*3,0,sigma)
#  f0 = dualpathSvd2(y0,D,verbose=F, maxsteps = 10)
#  # collect
#  for(maxsteps in 1:5)  print(head(sort(getGammat(f0,y0,5)%*%y0)))
#}

