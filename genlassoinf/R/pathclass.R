##' Constructor should take in the bare minimum to form an empty object
path <- function(...) {
    emptypath = list(lambda="", beta = "", fit="", Gammat="", u="", df="",
                     ss="", states="", completepath="", bls="", pathobjs="",
                     nk="", action="", tab="")
    structure(emptypath,
              class = "path")
}

##' Check if object is of class "path"
is.path <- function(someobj){ inherits(someobj, "path") }

##' Form generic
print <- function(x,...) UseMethod("print")

##' Plotting function for "sim" class
##' @import hexbin
print.path = function(mypath){
    ## Print the directions and signs
    cat("The selection sequence is: ", fill=TRUE)
    cat(mypath$action, fill=TRUE)

    ## Print the final variable set
    cat("The final variables selected after step ", mypath$pathobjs$k, " and lam= ", mypath$lambda, fill=TRUE)
    cat(getB.from.actions(mypath$action), fill=TRUE)

    ## Print the D matrix used 
    cat("A glimpse of the D matrix is", fill=TRUE)
    Dmat = mypath$D
    print(Dmat[1:(min(nrow(Dmat),5)), 1:(min(ncol(Dmat),5))])
} 

