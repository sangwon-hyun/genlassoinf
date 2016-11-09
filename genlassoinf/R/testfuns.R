##'  Function to run post-selection tests on graphs and store p-values results
##'  Currently only is able to do 3-clusters
##'  Input : path object, y0, initial graph, maximum steps to take, Dmat, cluster size, cluster #
##'  Output: 6 p values for each possible segment test.
sbmsim = function(f0,y0,mygraph, maxsteps,Dmat,clustersize,nclust=3,verbose=FALSE, stoptime=F){
  stopifnot(nclust==3)
  pvec = rep(NA,6)
  testtime.vec = rep(NA,6)
  
  # set initial graph's grouping
  prev.groups = list()
  for(member in unique(igraph::clusters(mygraph)$membership))  prev.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
  
  prev.nclust = length(prev.groups)

  # do clustering and conduct tests along the way (and plot them)
  for(step in 1:maxsteps){
  
    
    if(verbose){  cat(step , "out of", maxsteps, fill = TRUE) }

    # Get information about that step
      G0 = getGammat(f0, y0, step)
      Dmat.curr = Dmat[-getB.from.actions(f0$action[1:step]),]
      mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat.curr), mode="undirected")
      if(nrow(rbind(Dmat.curr))==0) break
      

    # Keep track of connected components
      curr.groups = list()
      vsegment = rep(0,length(y0))
      
      for(member in unique(igraph::clusters(mygraph)$membership)){ 
        curr.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
      }

    curr.nclust = length(curr.groups)
    # if cluster was created
    if(curr.nclust == prev.nclust + 1){
      if(verbose==TRUE) { print("cluster created!") }
      matched.curr.group.ind = c()
      # Which cluster was created? Scan all groups from prev and curr; what has changed?
      for(jj in 1:prev.nclust){
        for(kk in 1:curr.nclust){
          if( setequal(prev.groups[[jj]], curr.groups[[kk]])){
            vsegment[prev.groups[[jj]]] = 0 # if _any_ of the previous groups match with current groups exactly, assign zero
            matched.curr.group.ind = c(matched.curr.group.ind, kk)
          }
        }
      }
      # create segment test vector
      curr.test.group.ind = (if(length(matched.curr.group.ind) == 0) (1:curr.nclust) else (1:curr.nclust)[-matched.curr.group.ind])
      group1 = curr.groups[[curr.test.group.ind[1]]]
      group2 = curr.groups[[curr.test.group.ind[2]]]
      vsegment[group1] = -1/length(group1)
      vsegment[group2] =  1/length(group2)
                
      
      # conduct segment test
      pseg = pval.fl1d(y0, G0, vsegment, sigma, approxtype="rob")

      # collect p-value only if it corresponds to one of the pairs!
      pair.list = list(list(1:clustersize, (1:clustersize)+clustersize ),
                      list((1:clustersize)+clustersize,(1:clustersize)+2*clustersize ),
                      list(1:clustersize,(1:clustersize)+2*clustersize),
                      list(1:(2*clustersize),(1:clustersize)+2*clustersize),
                      list(c(1:clustersize,(1:clustersize)+2*clustersize),(1:clustersize)+clustersize),
                      list(1:clustersize,(1+clustersize):(3*clustersize)))
      pair.list = lapply(pair.list,function(lst){ list(as.integer(lst[[1]]), as.integer(lst[[2]] ))})
                     
                      
      # record the p-value in the place where the test is defined!
      for(ll in 1:length(pair.list)){
         if( all(list(group1, group2) %in% pair.list[[ll]])){
           pvec[ll] = pseg
           testtime.vec[ll] = step
         }
      }
      
    } else if (curr.nclust +1 == prev.nclust){
      # if clusters was merged
       warning("cluster being CREATED has not been coded yet")
    } else {
      if(verbose==TRUE) { print("no cluster created!") }
    }
           
    # update groups for next iteration
    prev.groups = curr.groups
    prev.nclust = curr.nclust
    
    # break if there are enough clusters
#      if(precurr.nclust > 4){
#        cat("found enough clusters, exiting loop")
#        break
#      }
  }
  
  if(stoptime==T){
    return(list(pvec=pvec,testtime.vec=testtime.vec))
  } else {
    return(pvec)
  }
}





##' Runs post-selection tests on 2d graphs and store p-values , and outputs 6 p
##' values for each possible segment test.
##' @param f0 path object
##' @param y0 data
##' @param mygraph initial graph
##' @param maxsteps maximum steps to take
##' @param Dmat graph penalty matrix to be used
##' @param clustersize size of cluster
##' @param verbose whether to be loud or not
graph2dsim <- function(f0,y0,mygraph, maxsteps,Dmat,clustersize,verbose=FALSE){
    results = list()
  
    ## set initial graph's grouping
    prev.groups = list()
    for(member in unique(igraph::clusters(mygraph)$membership)){
        prev.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
    }
  
    prev.nclust = length(prev.groups)
    count = 1
  
    # Do clustering and conduct tests along the way (and plot them)
    for(step in 1:maxsteps){
      if(verbose){  cat(step , "out of", maxsteps, fill = TRUE) }
      
      # Get information about that step
      G0 = getGammat(f0, y0, step)
      Dmat.curr = Dmat[-getB.from.actions(f0$action[1:step]),]
      mygraph = graph_from_adjacency_matrix(getadjmat.from.Dmat(Dmat.curr), mode="undirected")
      if(nrow(rbind(Dmat.curr))==0) break
      
  
      # Keep track of connected components
      curr.groups = list()
      vsegment = rep(0,length(y0))
      
      for(member in unique(igraph::clusters(mygraph)$membership)){ 
        curr.groups[[member]] = which(igraph::clusters(mygraph)$membership == member)
      }
      curr.nclust = length(curr.groups)
      
      # if cluster was created
      if(curr.nclust == prev.nclust + 1){
        if(verbose==TRUE) { print("cluster created!") }
        matched.curr.group.ind = c()
        # Which cluster was created? Scan all groups from prev and curr; what has changed?
        for(jj in 1:prev.nclust){
          for(kk in 1:curr.nclust){
            if( setequal(prev.groups[[jj]], curr.groups[[kk]])){
              vsegment[prev.groups[[jj]]] = 0 # if _any_ of the previous groups match with current groups exactly, assign zero
              matched.curr.group.ind = c(matched.curr.group.ind, kk)
            }
          }
        }
        # create segment test vector
        curr.test.group.ind = (if(length(matched.curr.group.ind) == 0) (1:curr.nclust) else (1:curr.nclust)[-matched.curr.group.ind])
        group1 = curr.groups[[curr.test.group.ind[1]]]
        group2 = curr.groups[[curr.test.group.ind[2]]]
        
        vsegment[group1] = +1/length(group1)
        vsegment[group2] =  -1/length(group2)    
        
        # adjust the sign of vsegment
        if(vsegment%*%y0 <= 0) vsegment = -vsegment 
        
        # conduct segment test
        pseg = pval.fl1d(y0, G0, vsegment, sigma, approxtype="rob")
  
        # record p-value and the groups that are being tested (and parse them later when aggregating results)
        results[[count]] = list(p = pseg, testttime = step, groups = list(group1=group1,
                                                            group2=group2))
        count = count + 1
        
      } else if (curr.nclust +1 == prev.nclust){
        # if clusters was merged
         warning("cluster being CREATED has not been coded yet")
      } else {
        if(verbose==TRUE) { print("no cluster created!") }
      }
             
      # update groups for next iteration
      prev.groups = curr.groups
      prev.nclust = curr.nclust
    }
    return(results)
}
