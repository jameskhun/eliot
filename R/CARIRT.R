#' A fast C++ implementation of the conditional autoregressive scaling model
#' 
#' 
#' \code{CARIRT} 
#' @importFrom Rcpp evalCpp
#' @useDynLib eliot
#' @param ConnectionMatrix
#' @param VoteMatrix 
#' @param IRT.polarity
#' @param burn.iter
#' @param mcmc.iter
#' @param thin.iter
#' @param store.alpha
#' @param store.beta
#' @param store.lambda
#' @param store.average.pref
#' @param store.latent.pref
#' @param store.expected.pref
#' @export
CARIRT <- function(ConnectionMatrix,VoteMatrix,IRT.polarity=1,burn.iter=100,mcmc.iter=100,thin.iter=1,store.alpha=FALSE,store.beta=FALSE,store.lambda=FALSE,store.average.pref=FALSE,store.latent.pref=FALSE,store.expected.pref=FALSE){
	
library(Rcpp)
library(RcppArmadillo)

if (dim(ConnectionMatrix)[1] != dim(ConnectionMatrix)[2]) stop("Matrix of Connections between Votes must be Square.")
if (dim(ConnectionMatrix)[1] != dim(VoteMatrix)[1]) stop("Matrix of Connections between Votes must have same dimensions as the Number of Votes.")
if (sum(ConnectionMatrix < 0) > 0) stop("Cannot have negative elements in matrix of connections between votes.")
if (dim(VoteMatrix)[2] < abs(IRT.polarity)) stop("Invalid polarity for IRT, selected voter index is greater than Vote Matrix dimension.")

nvotes <- dim(VoteMatrix)[1]
nvoters <- dim(VoteMatrix)[2]

vote.names <- rownames(VoteMatrix)
if (is.null(vote.names)) vote.names <- paste("Vote",1:nvotes)
voter.names <- colnames(VoteMatrix)
if (is.null(voter.names)) voter.names <- paste("Voter",1:nvoters)
	
savechains <- c(store.alpha,store.beta,store.lambda,store.average.pref,store.expected.pref,store.latent.pref)
		
library(Rcpp)
library(RcppArmadillo)	
library(inline)

# Get Starting Values for Ideal Points by Eigenvalue Decomposition

print("Getting Starting Values for Ideal Points by Eigenvalue Decomposition...")

x <- t(VoteMatrix)
row.mean <- apply(x, 1, mean, na.rm=TRUE)
col.mean <- apply(x, 2, mean, na.rm=TRUE)
dc1 <- sweep(x, 1, row.mean)
dc2 <- sweep(dc1, 2, col.mean)
dc <- dc2 + mean(x, na.rm = T)
r <- cor(t(dc),use="pairwise")
r[is.na(r)] <- 0
e <- eigen(r)
v <- e$vectors[,1]
start.average.prefs <- (v - mean(v))/sd(v)
  
# Set correct polarity for starting values

specified.voter <- abs(IRT.polarity)
if (sign(start.average.prefs[specified.voter]) != sign(IRT.polarity)) start.average.prefs <- -start.average.prefs

# Generate Starting Values for Vote Parameters

nonUnanimous <- rowMeans(VoteMatrix,na.rm=TRUE) < 1
nonUnanimous <- replace(nonUnanimous,is.na(nonUnanimous),-1)
VoteMatrix <- replace(VoteMatrix,VoteMatrix == 0,-1)
VoteMatrix <- replace(VoteMatrix,is.na(VoteMatrix),0)
start.alpha  <- rep(numeric(1),nvotes)
start.beta <- rep(integer(1),nvotes)
for (j in 1:nvotes){
	if (nonUnanimous[j] == 1){
		mean.maj <- mean(start.average.prefs[VoteMatrix[j,] == 1])
		mean.dis <- mean(start.average.prefs[VoteMatrix[j,] == -1])
		start.alpha[j] <- (mean.maj + mean.dis)/2
		start.beta[j] <- sign(mean.maj - mean.dis)
	} 
	if (nonUnanimous[j] == 0){
		if (runif(1) > 0.5){
			start.alpha[j] <- min(start.average.prefs[VoteMatrix[j,] == 1])
			start.beta[j] <- 1
		} else {
			start.alpha[j] <- max(start.average.prefs[VoteMatrix[j,] == 1])
			start.beta[j] <- -1
		}	
	}
	if (nonUnanimous[j] == -1){
		if (runif(1) > 0.5){
			start.alpha[j] <- rnorm(1)
			start.beta[j] <- 1
		} else {
			start.alpha[j] <- rnorm(1)
			start.beta[j] <- -1
		}	
	}
}

# Standardize alpha start values
start.alpha <- (start.alpha - mean(start.alpha))/sd(start.alpha)

# Prepare data for passing to C++ 

VoteMatrixdims <- dim(VoteMatrix)
VoteMatrix <- as.integer(VoteMatrix)
dim(VoteMatrix) <- VoteMatrixdims

Output <- CARIRTcpp(ConnectionMatrix,VoteMatrix, nonUnanimous, start.average.prefs, start.alpha, start.beta,savechains,burn.iter,mcmc.iter,thin.iter)

print("C++ Code Completed")

# Set appropriate dimensions for chain objects

dim(Output$alpha.chain) <- c(nvotes,mcmc.iter^store.alpha)
dimnames(Output$alpha.chain) <- c(list(vote.names,1:(mcmc.iter^store.alpha)))

dim(Output$beta.chain) <- c(nvotes,mcmc.iter^store.beta)
dimnames(Output$beta.chain) <- c(list(vote.names,1:(mcmc.iter^store.beta)))

dim(Output$lambda.chain) <- c(1,mcmc.iter^store.lambda)
dimnames(Output$lambda.chain) <- c(list(c("psi"),1:(mcmc.iter^store.lambda)))

dim(Output$average.pref.chain) <- c(nvoters,mcmc.iter^store.average.pref)
dimnames(Output$average.pref.chain) <- c(list(voter.names,1:(mcmc.iter^store.average.pref)))

dim(Output$expected.pref.chain) <- c(nvotes,nvoters,mcmc.iter^store.expected.pref)
dimnames(Output$expected.pref.chain) <- c(list(vote.names,voter.names,1:(mcmc.iter^store.expected.pref)))

dim(Output$latent.pref.chain) <- c(nvotes,nvoters,mcmc.iter^store.latent.pref)
dimnames(Output$latent.pref.chain) <- c(list(vote.names,voter.names,1:(mcmc.iter^store.latent.pref)))

dim(Output$ll.chain) <- mcmc.iter

class(Output) <- "CARIRT"

return(Output)
	
}