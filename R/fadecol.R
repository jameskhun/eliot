#' Function to fade colours by specified alpha
#' 
#' 
#' \code{fadecol} takes existing colours and applies user-specified alpha levels to them.
#' @param colours a vector providing any of the three kinds of R color specifications: hex, character string, or integer.
#' @param alpha the desired alpha level, must be in the interval 0 to 1. 
#' @export	
fadecol <- function(colours,alpha=1){
	rgb.dat <- col2rgb(colours)
	return(rgb(rgb.dat[1,],rgb.dat[2,],rgb.dat[3,],alpha*255,maxColorValue=255))
	}	