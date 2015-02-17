#' Convenience function to round values for display.
#' 
#' 
#' \code{round1}, \code{round2}, and \code{round2} round to 1, 2, and 3 points after the decimal, using the \code{formatC} function to retain trailing zeroes.
#' @param x 
#' @export
round1 <- function(x) formatC(x,digits=1,format="f")

#' Convenience function to round values for display.
#' 
#' 
#' \code{round1}, \code{round2}, and \code{round2} round to 1, 2, and 3 points after the decimal, using the \code{formatC} function to retain trailing zeroes.
#' @param x 
#' @export
round2 <- function(x) formatC(x,digits=2,format="f")

#' Convenience function to round values for display.
#' 
#' 
#' \code{round1}, \code{round2}, and \code{round2} round to 1, 2, and 3 points after the decimal, using the \code{formatC} function to retain trailing zeroes.
#' @param x 
#' @export
round3 <- function(x) formatC(x,digits=3,format="f")

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


#' Convenience functions for binary data transformations.
#' 
#' 
#' \code{logistic}.
#' 
#' @param x 
#' @export
logistic <- function(x) exp(x)/(1+exp(x))

#' Convenience functions for binary data transformations.
#' 
#' 
#' \code{logit}.
#' 
#' @param p 
#' @export
logit <- function(p) log((p)/(1-p))

#' Convenience functions for binary data transformations.
#' 
#' 
#' \code{invprobit}.
#' 
#' @param x 
#' @export
invprobit <- function(x) pnorm(x,0,1)

#' Convenience functions for binary data transformations.
#' 
#' 
#' \code{logit}.
#' 
#' @param p 
#' @export
probit <- function(p) qnorm(p,0,1)

#' Create empty plot.
#' 
#' 
#' \code{emptyplot}.
#' 
#' @param p 
#' @export
emptyplot <- function(xlim=c(-1,1),ylim=c(-1,1),xlab="",ylab="",main="",sub="",...) plot(0,0,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,sub=sub,axes=FALSE,...)

#' Add shaded confidence/prediction bands to a plot.
#' 
#' 
#' \code{shadeBands}.
#' 
#' @param x values of x at which interval has been calculated
#' @param lo lower bound of interval, at each value of x
#' @param hi upper bound of interval, at each value of x
#' @param ... additional arguments for the polygon function.
#' @export
shadeBands <- function(x,lo,hi,...) polygon(c(x,rev(x)),y=c(lo,rev(hi)),col=...,border=NA)

