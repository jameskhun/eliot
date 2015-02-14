#' A fast C++ implementation of the Wordfish text model
#' 
#' 
#' \code{wordminnow} implements Slapin and Proksch's Wordfish model in C++ for faster implementation.  Based on Will Lowe's EM algorithm R code in the package \code{austin}.
#' @importFrom Rcpp evalCpp
#' @useDynLib eliot
#' @param wfm The word/feature frequency matrix on which the model will be fit.
#' @param dir 
#' @param priors 
#' @param tol 
#' @export
wordminnow <- function(wfm, dir=c(1, 10), control=list(priors = c(Inf,Inf,3,1), tol = c(1e-6,1e-8))) {
	
  return(wordfishcpp(wfm, as.integer(dir), 1/(control$priors^2), control$tol))
    
}


