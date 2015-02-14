#' Function to standardize a vector
#' 
#' 
#' \code{standardize} subtracts the mean from and divides by the standard deviation of a vector.
#' @param x 
#' @export
standardize <- function(x,na.rm=FALSE) return((x-mean(x,na.rm=na.rm))/sd(x,na.rm=na.rm))
