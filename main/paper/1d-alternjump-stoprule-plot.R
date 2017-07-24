## Make sure you're working from [dropboxfolder]/code
source("settings.R")
source('helpers.R') # generic helpers
source('funs.R')
source('testfuns.R')
source('dualPathSvd2.R')
library(genlasso)
library(RColorBrewer)

outputdir="output"
codedir="."

################
### n = 20 #####
################

#

# Aggregate everything for n=20 (slower but readable)
# The function getpowers.from.chunk() is in testfuns.R, and reads simulation results
  powers.20 = list(powers.bic  = getpowers.from.chunks(20, "verdicts.bic"),
                   powers.bic2 = getpowers.from.chunks(20, "verdicts.bic2"),
                   powers.ebic = getpowers.from.chunks(20, "verdicts.ebic"),
                   powers.sbic = getpowers.from.chunks(20, "verdicts.sbic"),
                   powers.aic = getpowers.from.chunks(20, "verdicts.aic"),
                   powers.fixed1 = getpowers.from.chunks(20, "verdicts.fixed1"),
                   powers.fixed2 = getpowers.from.chunks(20, "verdicts.fixed2"),
                   powers.oracle = getpowers.oracle(20))
                   
  stoptimes.20 = getstoptimes(n=20,ngrain=20,nsim=100000,nchunks=100,file.id="bic-onejump-segmentsize-allbics")
#  n=20
#  ngrain=20
#  nsim = 100000
#  nchunks = 100
#  sigmalist = seq(from=0.1,to=3,length=20)
#  file.id="bic-onejump-segmentsize-allbics"
 # save(powers.20, stoptimes.20, ngrain,nsim,nchunks,sigmalist, file=file.path(outputdir, "bic-info-20.Rdata"))
  load(file=file.path(outputdir, "powers-20.Rdata"))

## make plot of powers
  pdf(file.path(outputdir,paste0( file.id,"-largesim", n/2, ".pdf")), width=10, height=6)
    makeplot(20, sigmalist, powers.20, stoptimes.20)
  dev.off()




