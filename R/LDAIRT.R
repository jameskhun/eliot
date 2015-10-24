#' A fast C++ implementation of the LDA + IRT scaling model
#' 
#' 
#' \code{LDAIRT} 
#' @importFrom Rcpp evalCpp
#' @useDynLib eliot
#' @param WordSparseMatrix
#' @param VoteMatrix 
#' @param LDA.topics
#' @param LDA.alpha
#' @param LDA.beta
#' @param IRT.discrim.prec
#' @param IRT.polarity
#' @param burn.iter
#' @param mcmc.iter
#' @param save.word.chain
#' @param save.topic.chain
#' @param save.ideal.chain
#' @param save.discrim.chain
#' @export
LDAIRT <- function(WordSparseMatrix,VoteMatrix,LDA.topics=10,LDA.alpha=1/LDA.topics,LDA.beta=1,IRT.discrim.prec=0.25,IRT.polarity=1,burn.iter=100,mcmc.iter=100,save.word.chain=FALSE,save.topic.chain=FALSE,save.ideal.chain=TRUE,save.discrim.chain=FALSE){
	
library(Rcpp)
library(RcppArmadillo)

if (WordSparseMatrix$nrow != dim(VoteMatrix)[1]) stop("Word Sparse Matrix and Vote Matrix must have the same number of rows.")
if (LDA.topics > WordSparseMatrix$nrow/5) stop("Too many topics requested, fewer than 5 votes per topic.")
if (dim(VoteMatrix)[2] < abs(IRT.polarity)) stop("Invalid polarity for IRT, selected voter index is greater than Vote Matrix dimension.")

ntopics <- LDA.topics
nterms <- WordSparseMatrix$ncol
nwords <- sum(WordSparseMatrix$v)
nvotes <- dim(VoteMatrix)[1]
nvoters <- dim(VoteMatrix)[2]
ndocterms <- length(WordSparseMatrix$v)

term.names <- WordSparseMatrix$dimnames$Terms
vote.names <- rownames(VoteMatrix)
if (is.null(vote.names)) vote.names <- paste("Vote",1:nvotes)
voter.names <- colnames(VoteMatrix)
if (is.null(voter.names)) voter.names <- paste("Voter",1:nvoters)
topic.names <- paste("Topic",1:ntopics)

SaveChains <- c(save.word.chain,save.topic.chain,save.ideal.chain,save.discrim.chain)		
	
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
r <- suppressWarnings(cor(t(dc),use="pairwise"))
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

WordMatrixDocs <- WordSparseMatrix$i - 1
WordMatrixTerms <- WordSparseMatrix$j - 1
WordMatrixCounts <- WordSparseMatrix$v

Output <- LDAIRTcpp(WordMatrixDocs, WordMatrixTerms, WordMatrixCounts,VoteMatrix,ntopics,nterms,nwords, ndocterms,LDA.alpha,LDA.beta,IRT.discrim.prec,start.ideal.points,c(save.word.chain, save.topic.chain, save.ideal.chain, save.discrim.chain),burn.iter,mcmc.iter)

print("C++ Code Completed")

# Set appropriate dimensions for chain objects

dim(Output$word.chain) <- c(ntopics,nterms,mcmc.iter^save.word.chain)
dimnames(Output$word.chain) <- c(list(topic.names,term.names,1:(mcmc.iter^save.word.chain)))
dim(Output$topic.chain) <- c(ntopics,nvotes,mcmc.iter^save.topic.chain)
dimnames(Output$topic.chain) <- c(list(topic.names,vote.names,1:(mcmc.iter^save.topic.chain)))
dim(Output$ideal.chain) <- c(ntopics,nvoters,mcmc.iter^save.ideal.chain)
dimnames(Output$ideal.chain) <- c(list(topic.names,voter.names,1:(mcmc.iter^save.ideal.chain)))
dim(Output$discrim.chain) <- c(nvotes,2,mcmc.iter^save.discrim.chain)
dimnames(Output$discrim.chain) <- c(list(vote.names,c("IRTalpha","IRTbeta"),1:(mcmc.iter^save.discrim.chain)))
dim(Output$rho.chain) <- mcmc.iter
dim(Output$ll.lda.chain) <- mcmc.iter
dim(Output$ll.irt.chain) <- mcmc.iter

dim(Output$dbar) <- 1
dim(Output$pd) <- 1
dim(Output$dic) <- 1

Output$dbar.error <- 2*sd(Output$ll.irt.chain)/sqrt(mcmc.iter)


print(paste("IRT Deviance Information Criterion:",round(Output$dic,2)))

gewekep.lda <- t.test(Output$ll.lda.chain[1:(trunc(mcmc.iter/4))],Output$ll.lda.chain[(trunc(3*mcmc.iter/4)):mcmc.iter])$p.value
gewekep.irt <- t.test(Output$ll.irt.chain[1:(trunc(mcmc.iter/4))],Output$ll.irt.chain[(trunc(3*mcmc.iter/4)):mcmc.iter])$p.value

print(paste("Geweke p-value for equal mean LDA log-likelihood in first and last quarter of simulation:",round(gewekep.lda,6)))
print(paste("Geweke p-value for equal mean IRT log-likelihood in first and last quarter of simulation:",round(gewekep.irt,6)))

class(Output) <- "LDAIRT"

return(Output)
	
}