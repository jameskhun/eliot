LIRT <- function(LambdaMatrix,VoteMatrix,IRT.discrim.prec=0.25,IRT.polarity=1,burn.iter=100,mcmc.iter=100,save.ideal.chain=TRUE,save.discrim.chain=FALSE){
	
library(Rcpp)
library(RcppArmadillo)

if (dim(LambdaMatrix)[1] != dim(VoteMatrix)[1]) stop("Lambda Matrix and Vote Matrix must have the same number of rows.")
if (dim(VoteMatrix)[2] < abs(IRT.polarity)) stop("Invalid polarity for IRT, selected voter index is greater than Vote Matrix dimension.")

ntopics <- dim(LambdaMatrix)[2]
nvotes <- dim(VoteMatrix)[1]
nvoters <- dim(VoteMatrix)[2]

vote.names <- rownames(VoteMatrix)
if (is.null(vote.names)) vote.names <- paste("Vote",1:nvotes)
voter.names <- colnames(VoteMatrix)
if (is.null(voter.names)) voter.names <- paste("Voter",1:nvoters)
topic.names <- colnames(LambdaMatrix)
if (is.null(topic.names)) topic.names <- paste("Topic",1:ntopics)

SaveChains <- c(TRUE,TRUE,save.ideal.chain,save.discrim.chain)		
	
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
start.ideal.points <- (v - mean(v))/sd(v)
  
# Set correct polarity for starting values

specified.voter <- abs(IRT.polarity)
if (sign(start.ideal.points[specified.voter]) != sign(IRT.polarity)) start.ideal.points <- -start.ideal.points

# Prepare data for passing to C++ 

VoteMatrixdims <- dim(VoteMatrix)
VoteMatrix <- as.integer(VoteMatrix)
dim(VoteMatrix) <- VoteMatrixdims
VoteMatrix <- replace(VoteMatrix,is.na(VoteMatrix),-9)


Output <- LIRTcpp(LambdaMatrix,VoteMatrix,IRT.discrim.prec,start.ideal.points,SaveChains,burn.iter,mcmc.iter)

print("C++ Code Completed")

# Set appropriate dimensions for chain objects

dim(Output$ideal.chain) <- c(ntopics,nvoters,mcmc.iter^save.ideal.chain)
dimnames(Output$ideal.chain) <- c(list(topic.names,voter.names,1:(mcmc.iter^save.ideal.chain)))
dim(Output$discrim.chain) <- c(nvotes,2,mcmc.iter^save.discrim.chain)
dimnames(Output$discrim.chain) <- c(list(vote.names,c("IRTalpha","IRTbeta"),1:(mcmc.iter^save.discrim.chain)))
dim(Output$rho.chain) <- mcmc.iter
dim(Output$ll.irt.chain) <- mcmc.iter
dim(Output$dbar) <- 1
dim(Output$pd) <- 1
dim(Output$dic) <- 1

Output$dbar.error <- 2*sd(Output$ll.irt.chain)/sqrt(mcmc.iter)

print(paste("IRT Deviance Information Criterion:",round(Output$dic,2)))

gewekep.irt <- t.test(Output$ll.irt.chain[1:(trunc(mcmc.iter/4))],Output$ll.irt.chain[(trunc(3*mcmc.iter/4)):mcmc.iter])$p.value

print(paste("Geweke p-value for equal mean IRT log-likelihood in first and last quarter of simulation:",round(gewekep.irt,6)))

class(Output) <- "LIRT"

return(Output)
	
}