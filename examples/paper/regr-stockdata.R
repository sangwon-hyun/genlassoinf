## MAKE sure you're working from [dropboxfolder]/code
source("settings.R")
source('selectinf/selectiveInference/R/funs.inf.R')
source('funs.R')
source('dualPathSvd2.R')
library(genlasso)
library(Matrix)
library(RColorBrewer)
load(file="examples/stockData.RData")
outputdir = "output"
codedir="."

# Format data and response
center = function(vec){
  as.numeric(scale(vec, center=TRUE, scale=FALSE))
}  
log.serial.transform = function(vec){
  ratios =  vec[2:length(vec)] / vec[1:(length(vec)-1)]
  logratios = log(ratios)
  return(logratios)
}  
X.uncentered = tricky_prices[,c(1,2,3)]#c(1,2,3,4,5)
## matplot(X.uncentered)
## matplot(tricky_prices[,30],type='o',pch=sapply(1:30,toString))
X = apply(X.uncentered,2,log.serial.transform)
X = apply(X,2,center)
TT = nrow(X)
J = ncol(X)

## This setting with any noise detects a spurious cut
beta0 = c(rep(c(-10,10,-10),each=TT/3),-10,-10, 
          rep(c(-10,10),each=TT/2),10,
          rep(10,each=TT))/10

X.augmented = do.call(cbind, lapply(1:ncol(X), function(icol) diag(X[,icol])))
mn = X.augmented %*% beta0

#########################
## Run and save path  ###
#########################
declutter.window = 10
do.declutter = T
sigma = 0.01/5
maxsteps = 20
seed = 4
set.seed(seed)

y = mn + rnorm(TT,0,sigma)
el.net.penalty = 0.001^2
X.tilde = rbind(X.augmented,diag(rep(sqrt(el.net.penalty), ncol(X.augmented)) ))
y.tilde = c(y,rep(0,ncol(X.augmented)))

D0 = makeDmat(TT,order=0)
D.orig = do.call(bdiag,lapply(1:J,function(jj) D0))
ginvX.tilde = ginv(X.tilde)
ginvX.orig  = ginv(X)
D.tilde = as.matrix(D.orig %*% ginvX.tilde)

# Obtain path
f0 = dualpathSvd2(y.tilde, D.tilde, maxsteps=maxsteps,verbose=T)
actions=f0$action[f0$action!=0]
states = get.states(f0$action)
bic = getbic.regression(y0.orig=y, f0=f0, sigma=sigma, maxsteps = maxsteps-1,
                       X.orig = X, ginvX.orig = ginvX.orig, D.orig = D.orig)
consec=2
stop.time = which.rise(bic,consec=consec) -1
stop.time = pmin(stop.time,length(y)-consec-1)         

plot(bic); text((1:length(f0$action)+1), rep(0.25,length(f0$action)), labels=f0$action)
abline(v=stop.time+1)


Gobj.new.with.stoptime = getGammat.with.stoprule(obj=f0, y=y.tilde, y0.orig=y,
                                                 condition.step = stop.time+consec,
                                                 stoprule = "bic", sigma=sigma,
                                                 consec=consec, maxsteps=maxsteps, type="regression",
                                                 X.orig=X, ginvX.orig = ginvX.orig,D.orig=D.orig)
G = Gobj.new.with.stoptime$Gammat
u = Gobj.new.with.stoptime$u; #u = rep(0,nrow(G))


## ## Load and harvest some things from path
y = y.tilde[1:TT]  
states = get.states(f0$action)
final.model.orig = states[[stop.time+1]]
final.model.orig.signs = f0$ss[[stop.time+1]]

## Conduct Inference
if(min(G%*%y.tilde-u)<=0){ print("polyhedron is off");   print((G%*%y.tilde)[which(G%*%y.tilde<u)])}
final.model.decluttered = declutter(final.model.orig,declutter.window)
vsegments = list()
pseg = c()
pseg.lrt=c()
for(jj in 1:length(final.model.decluttered)){print(jj)
  test.knot = final.model.decluttered[jj]
  adj.knot  = final.model.decluttered[-jj]
  test.knot.sign = final.model.orig.signs[which(final.model.orig == test.knot)]
  vsegment = make.v.regression(test.knot, adj.knot, test.knot.sign, f0, TT, J, X, X.tilde)
  vsegment.long = c(vsegment, rep(0,TT*J))
  pseg[jj] = pval.fl1d(y=y.tilde, G=G, u=u, dik=vsegment.long, sigma=sigma)
}

round(pseg,3)
  
names(pseg)  = final.model.decluttered
betalist.ind = lapply(1:5, function(ii) (TT*(ii-1)+1):(TT*ii))
membership   = sapply(final.model.decluttered,
                      function(thiscut) which(unlist(lapply(betalist.ind,
                                                            function(thisvar.ind) any(thiscut %in% thisvar.ind)))))
membership.before.declutter = sapply(final.model.orig,
                                     function(thiscut) which(unlist(lapply(betalist.ind,
                                                                           function(thisvar.ind) any(thiscut %in% thisvar.ind)))))

## Save these results
save(maxsteps,D.tilde,D.orig,ginvX.tilde,ginvX.orig,beta0,X.augmented,beta0,X,
     TT,J,X.uncentered,y, consec,y.tilde,X.tilde, f0, bic, sigma, stop.time, G,
     u, seed, maxsteps,declutter.window, pseg, final.model.orig,
     final.model.orig.signs,final.model.decluttered,
     file=file.path(outputdir,paste0("stockData-sim.Rdata")))
