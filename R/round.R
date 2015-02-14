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