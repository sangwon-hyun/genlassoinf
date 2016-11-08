Synopsis
====================
Software to implement post-selection inference for the generalized lasso.

Code Example
====================
In R code, an example of a 1d-fused lasso inference would be:
```
# Generate data and obtain path for 3-step 1d fused lasso
sigma = 1
n = 40
mn = rep(c(0,5), each = n/2)
y  = mn + rnorm(n, 0, sigma)
D  = makeDmat(n, ord = 1)
mypath = dualpathSvd2(y, D, maxsteps, approx = T)
stoptime = 3
Gobj    = getGammat.naive(obj = mypath, y = y, condition.step=stoptime,
			maxsteps = maxsteps, D=D)

# Make segment test contrast for the first detected knot
d = getdvec(mypath, y, test.step, stoptime, type = "segment")

# Conduct TG test at the first knot
pval = pval.fl1d(y0, G, d, sigma,u)	
print(paste("p-value for segment test at location", f0$pathobj$B[1], "is", round(pval)))
```
Motivation
====================
Repeat text from Synopsis.

Installation
====================
*Before you proceed* You (the user) should manually set the variables `outputdir` to the directory you want to save things to, and `codedir` to your local directory of *code* folder. 

<!---
API Reference
====================
Not written yet.
-->

Organization and Tests
====================
General
-------
*funs.R*: Contains functions that implement the method.

*testfuns.R*: Contains one-time functions that were written specifically for simulations. If you're a user, you probably don't need to look at this.

*settings.R*: Obsolete; used to be like a config file that contains the user's choice of `outputdir` and `codedir`.



`Examples' folder
--------
The experiments that produce the simulation results in the paper.

### 1d fused lasso
*1d-intro.R*: Produces an introductory plot of the adjusted and unadjusted post-selection test (Figure 1)

*1d-testintro.R*: Produces an introductory plot of the spike and segment test contrast (Figure 3)

*1d-onejump.R*: Runs simulations and produces QQ plots for the one-jump example (Figure 6, top row)

*1d-alternjump.R*: Runs simulations and produces QQ plots for the alternating jump example (Figure 6, bottom row)

*1d-onejump-stoprule-sim.R*: Runs simulations for the one-jump example using stopping rules. (Figure 7)

*1d-onejump-stoprule-plot.R*: Plots this. 

*1d-alternjump-stoprule-sim.R*: Runs simulations for the alternating jump example, using stopping rules. (not in paper)

*1d-alternjump-stoprule-plot.R*: Plots this.

*1d-smuce.R*: Runs and plots SMUCE inference for a one-jump example. (Figure 8)

*1d-props.R*: Runs simple experiment in one-jump signal to see proportion of detections. (Caption of Figure 6)

### Trend filtering (linear)
*tf-onejump-sim.R*: Runs simulations for a simple linear trendfilter signal with one jump. (Figure 9)

*tf-onejump-plot.R*: Plots this.

*tf-example.R*: Produces introductory plot for linear trend filtering (Figure 4)



### Graph fused lasso (2d and general)
*graph-sim.R*: Runs simulations for a 3-group (10 node each) example in several signal-to-noise ratio settings, over many simulations. (not in paper now)

*graph-plot.R*: Plots this.

*graph-example.R*: Runs and plots a 3-group (20 node each) graph example. (Figure 5)

*2d-sim.R*: Runs simulations for 2d graph example with a 10x10 greyscale image (Figure 10)

*2d-plot.R*: Plots this

### Regression
*regr-stockdata-sim.R*: Runs DJIA stock data regression example (Figure 11)

*regr-stockdata-plot.R*: Plots this.


### Misc
*stepsign-plot.R*: Produces step-sign plots for 1d fused lasso and linear trend filtering. (Figure 13, only shows 1d fused lasso)

*declutter-plot.R*: Produces plots before and after decluttering, for 1d fused lasso and linear trend filtering (Figure 12)

*ic-plot.R*: Produces the information criterion (IC) plots in Figure 2.

### CGH data analysis
*cgh-sim.R*: Runs and saves the path results for the CGH data example (Figure 14)

*cgh-plot.R*: Plotting code for the CGH data example (Figure 14)

*cgh-signpic.R*: Plotting code for a step-sign plot for the CGH data example, with and without clustering. (Figure in appendix)

*cgh-samplesplitting.R*: Implements the sample splitting approach using the CGH data, and produces the p-value table. (Figure 14)

Contributors
====================
Ryan Tibshirani, Max G'Sell, Sangwon Hyun (Justin)

License
====================
Not written yet.
