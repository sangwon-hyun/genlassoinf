##' Constructor should take in the bare minimum to form an empty object
genlassoinf <- function(...) {
    emptypath = list(lambda="", beta = "", fit="", Gammat="", u="", df="",
                     ss="", states="", completepath="", bls="", pathobjs="",
                     shits.list="", nk="", action="", tab="", Gobj.naive="",
                     Gobj.stoprule = "")
    structure(emptypath,
              class = "path")
}

##' Check if object is of class "path"
is.genlassoinf <- function(obj){ inherits(obj, "genlassoinf") }

##' Plotting function for "path" class
##' @export
print.genlassoinf <- function(mypath){

    cat("The changepoints are", mypath$cp * mypath$cp.sign, fill=TRUE)

    ## ## Print the directions and signs
    ## cat("The selection sequence is: ", fill=TRUE)
    ## cat(mypath$action, fill=TRUE)

    ## ## Print the final variable set
    ## cat("The final variables selected after step ", mypath$pathobjs$k, "are", fill=TRUE)
    ## mytable = data.frame(matrix(nrow=1,ncol=length(mypath$lambda)+1))
    ## colnames(mytable) = c("lambda=", round(mypath$lambda,2))
    ## rownames(mytable) = ""
    ## mytable[1,2:ncol(mytable)] = getB.from.actions(mypath$action)
    ## mytable[1,1] = "knots="
    ## print(mytable)

    ## ## Print the D matrix used
    ## cat("A glimpse of the D matrix is", fill=TRUE)
    ## Dmat = mypath$D
    ## print(Dmat[1:(min(nrow(Dmat),5)), 1:(min(ncol(Dmat),5))])
}

##' Form generic
step_sign_plot <- function(x,...) UseMethod("step_sign_plot")

##' Takes a generalized lasso path object (from running \code{dualpathSvd2()})
##' and returns a numeric vector of values that you can plot as a step sign
##' plot, and optionally produces a step sign plot. 
##' @param obj product from running 1d fused lasso \code{dualpathSvd2()}.
##' @param stoptime Desired stop time of algorithm. Set to 0 for /no/
##'     changepoints, set to 1 for one changepoint, and so on.
##' @param postprocess TRUE if you want post-processing. Defaults to FALSE
##' @examples step.sign.plot(signs = c(+1,-1), final.model = c(20,36), n=60)
##' @export
step_sign_plot.genlassoinf <- function(obj, stoptime, postprocess=FALSE, plot=TRUE){

    ## Basic checks
    stopifnot(is.genlassoinf(obj))
    stopifnot(stoptime>0)

    ## Harvest things from obj
    signs = obj$ss[[stoptime+1]]
    final.model = obj$states[[stoptime+1]]
    s0 = step_sign_plot_inner(signs, final.model, n)

    ## Optionally, produce plot
    if(plot == TRUE){
        plot(s0, type='l', axes=F, xlab="", ylab="", ylim=range(s0) + c(-2, 2), lwd=2)
        axis(1)
        sign.chars = sapply(signs, function(mysign){if(mysign==+1)"+" else "-"})
        text(x=final.model, y=-0.15, label=sign.chars, cex=1.5)
        title(main="Step-sign plot")
        abline(v=final.model, col='lightgrey', lty=2)
    }
    return(s0)
}


##' Takes a generalized lasso path object (from running \code{dualpathSvd2()}) and
##' plots the mean and changepoints of a `good' model -- good is either the
##' model size selected by cross-validation, or selected by an IC rule.
##' @param obj product from running 1d fused lasso \code{dualpathSvd2()}.
##' @param stoprule Desired stopping rule of algorithm.
##' @export
plot.genlassoinf <- function(obj, stoprule = c("cv", "ic")){

    ## Basic things
    stoprule = match.arg(stoprule)
    if(stoprule == "ic"){
        if(!("stoprule" %in% objects(obj))){
           stop("Run stoptime() on your path object before trying to plot!")
        }
    }

    ## Plot settings (fixed for now)
    lcol.beta = "red"
    lty.beta = 1
    lwd.beta = 2
    lwd.abline = 2
    lty.abline = 2
    lcol.abline = 'lightgrey'
    xlab = "Coordinate"
    ylab = "y"
    pcol.dat = 'grey50'
    pch.dat = 16

    ## Make plot and add lines
    plot(y, axes=FALSE, xlab = xlab, ylab = ylab, col = pcol.dat, pch = pch.dat)
    axis(1);axis(2)
    abline(v=locs, col = lcol.abline, lty=lty.abline, lwd=lwd.abline)
    lines(obj$beta[,obj$stoptime], col = lcol.beta, lwd=lwd.fl)

    ## Add legend
    legend("topleft", legend = c("data", "primal fit", "loc"),
           lty = c(NA, lty.beta, lty.abline),
           col = c(pcol.dat, lcol.beta, lcol.abline), pch =c(pch.dat, NA,NA))
}
