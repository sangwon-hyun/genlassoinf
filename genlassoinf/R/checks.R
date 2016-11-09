##' Seeing if two Gamma matrices are the same, by matching them case-by-case
##' |Gamma1|, |Gamma2| are two Gamma matrices (not necessarily same number of rows)
##' y0 is the original response data, just used as a reference
checksame <- function(y0,Gamma1,Gamma2,noise=1){
  y0noisy   <- y0 + rnorm(length(y0), 0, noise) 
  return(all( Gamma1 %*% y0noisy >= 0 ) == all( Gamma2 %*% y0noisy >= 0 ))
} 

# Check correctness of polyhedron:
polyhedron.checks.out <- function(y0,G,u, throw.error=F){
  all.nonneg = all((G%*%y0 >= u))
  if(all.nonneg){ 
    return(all.nonneg) 
  } else {
    if(throw.error){ 
      stop("failed Gy >= u test")
    } else { 
      print("failed Gy >= u test")
    }
    return(all.nonneg)
  }
}        
