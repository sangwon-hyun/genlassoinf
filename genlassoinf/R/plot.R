##' Takes in the signs |signs| and the locations |final.model|
##'
##' @export
step.sign.plot <- function(signs, final.model, n){
    signs = signs[order(final.model)]
    signs = c(-signs[1], signs)
    cumul.signs = cumsum(signs)
    final.model = sort(final.model)
    nn = length(final.model)
    indices <- vector(mode = "list", length = nn+1)
      indices[2:nn] = Map(function(a,b){a:b},
                              final.model[1:(length(final.model)-1)]+1, 
                              final.model[2:length(final.model)])
      indices[[1]] = 1:final.model[1]
      indices[[nn+1]] = (final.model[length(final.model)]+1):n
    sign0 = do.call(c, 
             lapply(1:length(cumul.signs), 
                function(ii){rep(cumul.signs[ii], length(indices[[ii]]))})
                )
    return(sign0)
}

